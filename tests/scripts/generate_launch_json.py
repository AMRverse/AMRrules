#!/usr/bin/env python3
"""
generate_launch_json.py

Regenerates .vscode/launch.json from tests/example_commands.yaml, so your
debug configs can never silently drift from the CI/regression case list
again. Run this whenever you edit the manifest and want your VSCode debug
dropdown to stay in sync.

Debug configs write to tests/data/output (not example_output or
local_output) so debugging a single case never touches your baseline or
your regression-check scratch dir.

USAGE
  python tests/scripts/generate_launch_json.py
  python tests/scripts/generate_launch_json.py --include-debug-cases
"""

import argparse
import json
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
MANIFEST_PATH = REPO_ROOT / "tests" / "scripts" / "example_commands.yaml"
LAUNCH_JSON_PATH = REPO_ROOT / ".vscode" / "launch.json"
DEBUG_OUTPUT_DIR = "tests/data/local_output"


def case_to_config(case):
    args = [
        "--input", case["input"],
        "--output-prefix", case["name"],
        "--output-dir", DEBUG_OUTPUT_DIR,
    ]
    if case.get("organism"):
        args += ["--organism", case["organism"]]
    if case.get("organism_file"):
        args += ["--organism-file", case["organism_file"]]
    if case.get("sample_id"):
        args += ["--sample-id", case["sample_id"]]
    if case.get("flag_core"):
        args.append("--flag-core")
    if case.get("nr"):
        args += ["-nr", case["nr"]]
    args += case.get("extra_args", [])

    return {
        "name": case["name"],
        "type": "python",
        "request": "launch",
        "module": "amrrules",
        "cwd": "${workspaceFolder}",
        "console": "integratedTerminal",
        "args": args,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--include-debug-cases", action="store_true",
                         help="Also generate configs for cases under debug_cases in the manifest")
    args = parser.parse_args()

    manifest = yaml.safe_load(MANIFEST_PATH.read_text())
    cases = list(manifest.get("cases", []))
    if args.include_debug_cases:
        cases += manifest.get("debug_cases", [])

    launch_json = {
        "version": "0.2.0",
        "configurations": [case_to_config(c) for c in cases],
    }

    LAUNCH_JSON_PATH.parent.mkdir(parents=True, exist_ok=True)
    LAUNCH_JSON_PATH.write_text(json.dumps(launch_json, indent=4) + "\n")
    print(f"Wrote {len(cases)} debug configuration(s) to {LAUNCH_JSON_PATH}")


if __name__ == "__main__":
    main()
