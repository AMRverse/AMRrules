#!/usr/bin/env python3
"""
compare_outputs.py

Compare two (or many pairs of) TSV output files that share the same header,
to see exactly what changed between one version of your pipeline/code and
another. Rows are matched by a "key" (a combination of columns that
uniquely identifies a row), NOT by row position/order - so reordered rows
won't show up as spurious differences.

For each pair of files you get:
  - rows only in file A (dropped)
  - rows only in file B (added)
  - for rows present in both: which columns changed, and old -> new values

USAGE

  Single pair:
    python compare_outputs.py fileA.tsv fileB.tsv

  Single pair, custom key columns and output location:
    python compare_outputs.py fileA.tsv fileB.tsv \
        --key sample,drug,"drug class" \
        --out diff_report.tsv

  Batch mode - compare every file in dirA against the file of the same
  name in dirB (e.g. old_run/ vs new_run/), writing one report per file
  plus a summary table:
    python compare_outputs.py --batch old_run/ new_run/ --out-dir diffs/

If you don't pass --key, the script guesses a sensible key automatically
(see `guess_key_columns` below) and tells you what it picked.
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

# Columns that typically identify "what row this is" as opposed to columns
# that hold the actual result/value being tested. Different output formats
# from the pipeline use different column names, so we keep a list of
# "recipes" and use whichever one (a) exists in the file's columns and
# (b) actually uniquely identifies every row. Add new recipes here as new
# output formats come up - order matters, first match wins.
KEY_RECIPES = [
    # "genome_summary" style files: one row per sample/drug/drug-class
    ["sample", "drug", "drug class"],
    # "interpreted" style files: one row per genomic feature, exploded per
    # associated drug (so position alone isn't unique)
    ["Name", "variation type", "gene", "mutation", "ruleID", "drug", "drug class"],
]


def load_tsv(path):
    """Load a TSV as strings (so '0' vs 0.0 etc never causes false diffs),
    with NaNs turned into empty strings for clean comparison/printing."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return df


def guess_key_columns(df):
    """Try each recipe in KEY_RECIPES in order; return the first one whose
    columns are all present in df AND which uniquely identifies every row.
    Falls back to all columns if nothing matches."""
    for recipe in KEY_RECIPES:
        if all(c in df.columns for c in recipe) and not df.duplicated(subset=recipe).any():
            return recipe
    return list(df.columns)


def diff_two_files(path_a, path_b, key_cols=None, set_cols=None, ignore_cols=None):
    """Compare two TSVs with identical headers.

    ignore_cols: columns to exclude from comparison entirely (e.g. a
    'version' column that's expected to differ on every row and would
    otherwise drown out the diffs you actually care about).

    Returns a dict with:
      key_cols       - the key columns actually used
      added_rows     - DataFrame of rows only in B
      removed_rows   - DataFrame of rows only in A
      changed        - long-format DataFrame: key cols + column + old + new
      n_common_rows  - rows present in both, matched by key
      n_unchanged    - of those, how many are identical
    """
    df_a = load_tsv(path_a)
    df_b = load_tsv(path_b)

    if list(df_a.columns) != list(df_b.columns):
        raise ValueError(
            f"Column mismatch between {path_a} and {path_b}:\n"
            f"  A: {list(df_a.columns)}\n  B: {list(df_b.columns)}"
        )

    if key_cols is None:
        key_cols = guess_key_columns(df_a)

    for col in key_cols:
        if col not in df_a.columns:
            raise ValueError(f"Key column '{col}' not found in {path_a}")

    if df_a.duplicated(subset=key_cols).any() or df_b.duplicated(subset=key_cols).any():
        raise ValueError(
            f"Key columns {key_cols} do not uniquely identify rows in "
            f"{path_a} and/or {path_b}. Pass a more specific --key."
        )

    ignore_cols = set(ignore_cols or [])
    value_cols = [c for c in df_a.columns if c not in key_cols and c not in ignore_cols]

    merged = df_a.merge(
        df_b, on=key_cols, how="outer", suffixes=("_A", "_B"), indicator=True
    )

    added_rows = merged[merged["_merge"] == "right_only"][key_cols]
    removed_rows = merged[merged["_merge"] == "left_only"][key_cols]
    common = merged[merged["_merge"] == "both"]

    set_cols = set(set_cols or [])
    change_records = []
    for col in value_cols:
        col_a, col_b = f"{col}_A", f"{col}_B"
        sub = common[key_cols + [col_a, col_b]]

        if col in set_cols:
            def _norm(v):
                return frozenset(x.strip() for x in v.split(";") if x.strip())
            is_diff = sub[col_a].apply(_norm) != sub[col_b].apply(_norm)
        else:
            is_diff = sub[col_a] != sub[col_b]

        diffs = sub[is_diff]
        for _, row in diffs.iterrows():
            rec = {k: row[k] for k in key_cols}
            rec["column"] = col
            rec["old_value"] = row[col_a]
            rec["new_value"] = row[col_b]
            change_records.append(rec)

    changed = pd.DataFrame(change_records, columns=key_cols + ["column", "old_value", "new_value"])

    n_common_rows = len(common)
    n_changed_row_keys = changed[key_cols].drop_duplicates().shape[0] if not changed.empty else 0
    n_unchanged = n_common_rows - n_changed_row_keys

    return {
        "key_cols": key_cols,
        "added_rows": added_rows.reset_index(drop=True),
        "removed_rows": removed_rows.reset_index(drop=True),
        "changed": changed,
        "n_common_rows": n_common_rows,
        "n_unchanged": n_unchanged,
    }


def print_summary(label, result):
    print(f"\n=== {label} ===")
    print(f"Key columns used: {result['key_cols']}")
    print(f"Rows only in A (removed): {len(result['removed_rows'])}")
    print(f"Rows only in B (added):   {len(result['added_rows'])}")
    print(f"Rows in both:             {result['n_common_rows']} "
          f"({result['n_unchanged']} unchanged, "
          f"{result['n_common_rows'] - result['n_unchanged']} with at least one changed value)")
    if not result["changed"].empty:
        print("Changed value counts by column:")
        print(result["changed"]["column"].value_counts().to_string())


def write_report(result, out_path):
    """Write one combined, easy-to-scan TSV: a row per change, plus rows for
    added/removed keys, all tagged with a 'change_type' column."""
    key_cols = result["key_cols"]
    parts = []

    if not result["removed_rows"].empty:
        r = result["removed_rows"].copy()
        r["change_type"] = "row_removed"
        r["column"] = ""
        r["old_value"] = ""
        r["new_value"] = ""
        parts.append(r)

    if not result["added_rows"].empty:
        a = result["added_rows"].copy()
        a["change_type"] = "row_added"
        a["column"] = ""
        a["old_value"] = ""
        a["new_value"] = ""
        parts.append(a)

    if not result["changed"].empty:
        c = result["changed"].copy()
        c["change_type"] = "value_changed"
        parts.append(c)

    cols = key_cols + ["change_type", "column", "old_value", "new_value"]
    if parts:
        report = pd.concat(parts, ignore_index=True)[cols]
    else:
        report = pd.DataFrame(columns=cols)

    report.to_csv(out_path, sep="\t", index=False)
    return report


def run_single(path_a, path_b, key, out, set_cols, ignore_cols):
    result = diff_two_files(path_a, path_b, key_cols=key, set_cols=set_cols, ignore_cols=ignore_cols)
    print_summary(f"{Path(path_a).name}  vs  {Path(path_b).name}", result)
    report = write_report(result, out)
    print(f"\nFull diff written to: {out}  ({len(report)} rows)")
    return result


def run_batch(dir_a, dir_b, key, out_dir, pattern, set_cols, ignore_cols):
    dir_a, dir_b, out_dir = Path(dir_a), Path(dir_b), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files_a = {p.name: p for p in dir_a.glob(pattern)}
    files_b = {p.name: p for p in dir_b.glob(pattern)}
    common_names = sorted(set(files_a) & set(files_b))
    only_a = sorted(set(files_a) - set(files_b))
    only_b = sorted(set(files_b) - set(files_a))

    if only_a:
        print(f"Files only in {dir_a} (skipped): {only_a}")
    if only_b:
        print(f"Files only in {dir_b} (skipped): {only_b}")
    if not common_names:
        print("No matching filenames found in both directories.")
        return

    summary_rows = []
    for name in common_names:
        try:
            result = diff_two_files(files_a[name], files_b[name], key_cols=key, set_cols=set_cols, ignore_cols=ignore_cols)
        except ValueError as e:
            print(f"\n[SKIPPED] {name}: {e}")
            continue

        print_summary(name, result)
        out_path = out_dir / f"diff__{Path(name).stem}.tsv"
        write_report(result, out_path)

        summary_rows.append({
            "file": name,
            "rows_removed": len(result["removed_rows"]),
            "rows_added": len(result["added_rows"]),
            "rows_with_changes": result["n_common_rows"] - result["n_unchanged"],
            "rows_unchanged": result["n_unchanged"],
            "changed_columns": ", ".join(sorted(result["changed"]["column"].unique())) if not result["changed"].empty else "",
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = out_dir / "summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print(f"\n=== Batch summary written to {summary_path} ===")
    print(summary_df.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("path_a", nargs="?", help="First file (or, with --batch, first directory)")
    parser.add_argument("path_b", nargs="?", help="Second file (or, with --batch, second directory)")
    parser.add_argument("--key", help="Comma-separated key column(s), e.g. sample,drug,'drug class'. "
                                       "If omitted, guessed automatically.")
    parser.add_argument("--set-cols", help="Comma-separated columns where ';'-separated values should be "
                                            "compared as an unordered set (ignores order/exact duplicates), "
                                            "e.g. ruleIDs,'markers (non-S)'")
    parser.add_argument("--ignore-cols", help="Comma-separated columns to exclude from comparison entirely, "
                                               "e.g. version (useful for columns expected to differ on every "
                                               "row, like a tool/pipeline version stamp)")
    parser.add_argument("--out", default="diff_report.tsv", help="Output path for single-pair mode (default: diff_report.tsv)")
    parser.add_argument("--batch", action="store_true", help="Batch mode: path_a/path_b are directories")
    parser.add_argument("--out-dir", default="diffs", help="Output directory for batch mode (default: diffs/)")
    parser.add_argument("--pattern", default="*.tsv", help="Filename glob pattern for batch mode (default: *.tsv)")
    args = parser.parse_args()

    if not args.path_a or not args.path_b:
        parser.print_help()
        sys.exit(1)

    key = [c.strip().strip("'\"") for c in args.key.split(",")] if args.key else None
    set_cols = [c.strip().strip("'\"") for c in args.set_cols.split(",")] if args.set_cols else None
    ignore_cols = [c.strip().strip("'\"") for c in args.ignore_cols.split(",")] if args.ignore_cols else None

    if args.batch:
        run_batch(args.path_a, args.path_b, key, args.out_dir, args.pattern, set_cols, ignore_cols)
    else:
        run_single(args.path_a, args.path_b, key, args.out, set_cols, ignore_cols)


if __name__ == "__main__":
    main()
