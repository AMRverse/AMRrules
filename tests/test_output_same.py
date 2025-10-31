import sys
import json
import subprocess
import filecmp
from pathlib import Path
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
LAUNCH_JSON = PROJECT_ROOT / ".vscode" / "launch.json"
GOLDEN_DIR = PROJECT_ROOT / "tests" / "data" / "correct_output"
OUTPUT_DIR = PROJECT_ROOT / "tests" / "data" / "output"

# load launch.json at import time and build params
with open(LAUNCH_JSON, "r") as f:
    launch = json.load(f)

configs = launch.get("configurations", [])
params = []
ids = []
for cfg in configs:
    name = cfg.get("name", "unnamed")
    args = list(cfg.get("args", []))
    params.append((name, args))
    ids.append(name)

@pytest.mark.integration
@pytest.mark.parametrize("cfg_name,cfg_args", params, ids=ids)
def test_launch_config_runs(cfg_name, cfg_args, tmp_path):
    """
    Run each launch.json configuration against the CLI.
    Ensures the process exits 0 and expected outputs are created.
    If golden files exist in tests/data/correct_output they are compared.
    """
    args = list(cfg_args)  # work on a copy

    # make input path absolute if present
    if "--input" in args:
        i = args.index("--input")
        input_path = Path(args[i + 1])
        if not input_path.is_absolute():
            candidate = PROJECT_ROOT / input_path
            if candidate.exists():
                args[i + 1] = str(candidate)

    # ensure there is an output_prefix; if not, create one from config name
    if "--output_prefix" in args:
        pidx = args.index("--output_prefix")
        prefix = args[pidx + 1]
    else:
        prefix = cfg_name.replace(" ", "_")
        args += ["--output_prefix", prefix]

    # run the CLI with the same Python running pytest
    cmd = [sys.executable, "-m", "amrrules.cli"] + args
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise AssertionError(
            f"CLI failed for config '{cfg_name}' (rc={proc.returncode}):\n"
            f"CMD: {' '.join(cmd)}\n\nSTDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}"
        )

    # expected output file names used in your other tests
    interpreted = OUTPUT_DIR / f"{prefix}_interpreted.tsv"
    summary = OUTPUT_DIR / f"{prefix}_genome_summary.tsv"

    assert interpreted.exists(), f"Expected interpreted output missing: {interpreted}"
    assert summary.exists(), f"Expected summary output missing: {summary}"

    # if golden files exist, compare them byte-for-byte
    golden_interpreted = GOLDEN_DIR / interpreted.name
    golden_summary = GOLDEN_DIR / summary.name

    if golden_interpreted.exists():
        assert filecmp.cmp(str(golden_interpreted), str(interpreted), shallow=False), (
            f"Interpreted output differs for {cfg_name}"
        )
    if golden_summary.exists():
        assert filecmp.cmp(str(golden_summary), str(summary), shallow=False), (
            f"Summary output differs for {cfg_name}"
        )