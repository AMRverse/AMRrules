#!/usr/bin/env python3
"""
run_examples.py

Runs every case defined in tests/example_commands.yaml through amrrules,
writing results to a given --output-dir. This is the single execution path
used by:

  - CI (generate-example-outputs.yml), writing to tests/data/example_output
    -> this becomes the committed baseline "what output should look like"
  - your local regression check (tests/scripts/check_regressions.py), writing to
    a scratch dir, run against your current (possibly modified) code

USAGE

  # regenerate the committed baseline (what CI does)
  python tests/scripts/run_examples.py --output-dir tests/data/example_output

  # regenerate into a scratch dir using your current code, for comparison
  python tests/scripts/run_examples.py --output-dir tests/data/local_output

  # only run one or two cases while iterating on a fix
  python tests/scripts/run_examples.py --output-dir tests/data/local_output \
      --only test_ngono_20strains,test_kpneumo_MDR_nwtR

  # also run the debug-only cases (large/ungitted inputs - see manifest)
  python tests/scripts/run_examples.py --output-dir tests/data/local_output \
      --include-debug-cases
"""

import argparse
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST_PATH = REPO_ROOT / "tests" / "example_commands.yaml"


def load_cases(include_debug_cases=False):
    manifest = yaml.safe_load(MANIFEST_PATH.read_text())
    cases = list(manifest.get("cases", []))
    if include_debug_cases:
        cases += manifest.get("debug_cases", [])
    return cases


def build_command(case, output_dir):
    if bool(case.get("organism")) == bool(case.get("organism_file")):
        raise ValueError(
            f"Case '{case['name']}': specify exactly one of 'organism' or 'organism_file'"
        )

    cmd = [
        "amrrules",
        "--input", case["input"],
        "--output-prefix", case["name"],
        "--output-dir", str(output_dir),
    ]
    if case.get("organism"):
        cmd += ["--organism", case["organism"]]
    if case.get("organism_file"):
        cmd += ["--organism-file", case["organism_file"]]
    if case.get("sample_id"):
        cmd += ["--sample-id", case["sample_id"]]
    if case.get("flag_core"):
        cmd.append("--flag-core")
    if case.get("nr"):
        cmd += ["-nr", case["nr"]]
    cmd += case.get("extra_args", [])
    return cmd


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", required=True, help="Directory to write amrrules output to")
    parser.add_argument("--only", help="Comma-separated case names to run (default: all)")
    parser.add_argument("--include-debug-cases", action="store_true",
                         help="Also run cases listed under debug_cases in the manifest")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without running them")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cases = load_cases(include_debug_cases=args.include_debug_cases)
    if args.only:
        wanted = set(n.strip() for n in args.only.split(","))
        cases = [c for c in cases if c["name"] in wanted]
        missing = wanted - {c["name"] for c in cases}
        if missing:
            sys.exit(f"Unknown case name(s) in --only: {sorted(missing)}")

    print(f"Running {len(cases)} case(s) -> {output_dir}\n")

    failures = []
    for i, case in enumerate(cases, 1):
        cmd = build_command(case, output_dir)
        print(f"[{i}/{len(cases)}] {case['name']}")
        if args.dry_run:
            print("  " + " ".join(cmd))
            continue
        result = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
        if result.returncode != 0:
            failures.append(case["name"])
            print(f"  FAILED (exit {result.returncode})")
            print("  " + " ".join(cmd))
            if result.stderr.strip():
                print("  stderr (last 20 lines):")
                for line in result.stderr.strip().splitlines()[-20:]:
                    print("   ", line)

    if failures:
        print(f"\n{len(failures)} case(s) failed to run: {failures}")
        sys.exit(1)
    print("\nAll cases completed successfully.")


if __name__ == "__main__":
    main()
