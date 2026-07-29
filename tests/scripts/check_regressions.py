#!/usr/bin/env python3
"""
check_regressions.py

The "did I break anything?" one-liner. Regenerates example outputs using
your CURRENT (possibly modified) code, then diffs every file against the
committed baseline in tests/data/example_output using compare_outputs.py.

Before running this, make sure tests/data/example_output reflects main
(i.e. you haven't edited it locally) - e.g.:
    git fetch origin && git checkout origin/main -- tests/data/example_output

USAGE

  python tests/scripts/check_regressions.py

  # only check a subset of cases while iterating on a specific fix
  python tests/scripts/check_regressions.py --only test_ngono_20strains,test_kpneumo_MDR_nwtR

  # skip regenerating and just re-diff what's already in tests/data/local_output
  python tests/scripts/check_regressions.py --skip-generate

Output: one diff__<file>.tsv per output file under tests/data/regression_diffs/,
plus a summary.tsv you can scan in one go. Nothing here modifies
tests/data/example_output - that's your read-only baseline.
"""

import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
BASELINE_DIR = REPO_ROOT / "tests" / "data" / "example_output"
LOCAL_DIR = REPO_ROOT / "tests" / "data" / "local_output"
DIFF_DIR = REPO_ROOT / "tests" / "data" / "regression_diffs"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--only", help="Comma-separated case names to check (default: all)")
    parser.add_argument("--include-debug-cases", action="store_true")
    parser.add_argument("--skip-generate", action="store_true",
                         help="Don't regenerate; just diff tests/data/local_output as it currently stands")
    parser.add_argument("--ignore-cols", default="version",
                         help="Comma-separated columns to exclude from comparison (default: version)")
    args = parser.parse_args()

    if not BASELINE_DIR.exists():
        sys.exit(
            f"Baseline dir {BASELINE_DIR} doesn't exist. Pull it from main first, e.g.:\n"
            f"  git fetch origin && git checkout origin/main -- tests/data/example_output"
        )

    if not args.skip_generate:
        gen_cmd = [
            sys.executable, str(REPO_ROOT / "tests" / "scripts" / "run_examples.py"),
            "--output-dir", str(LOCAL_DIR),
        ]
        if args.only:
            gen_cmd += ["--only", args.only]
        if args.include_debug_cases:
            gen_cmd.append("--include-debug-cases")
        print("=== Regenerating outputs from current code ===")
        result = subprocess.run(gen_cmd, cwd=REPO_ROOT)
        if result.returncode != 0:
            sys.exit("Generation failed - fix the errors above before checking for regressions.")

    print("\n=== Comparing against baseline (tests/data/example_output) ===")
    diff_cmd = [
        sys.executable, str(REPO_ROOT / "tests" / "scripts" / "compare_outputs.py"),
        "--batch", str(BASELINE_DIR), str(LOCAL_DIR),
        "--out-dir", str(DIFF_DIR),
        "--ignore-cols", args.ignore_cols,
    ]
    subprocess.run(diff_cmd, cwd=REPO_ROOT)

    print(f"\nDone. Per-file diffs and summary.tsv are in {DIFF_DIR}")
    print("An empty diff__*.tsv (or a summary.tsv row with all zeros) means that file is unchanged.")


if __name__ == "__main__":
    main()
