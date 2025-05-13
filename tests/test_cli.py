import subprocess
import os

def test_cli():
    # Define paths
    input_file = "tests/data/input/test_ecoli_genome.tsv"
    rules_file = "tests/data/input/Ecoli_testRules_withCipro.txt"
    output_dir = "tests/data/output"
    output_prefix = "test_ecoli_genome"

    # Run the CLI
    result = subprocess.run(
        [
            "python", "src/amrrules/cli.py",
            "--input", input_file,
            "--output_dir", output_dir,
            "--output_prefix", output_prefix,
            "--rules", rules_file,
            "--annot_opts", "minimal"
        ],
        capture_output=True,
        text=True
    )

    # Assert the CLI ran successfully
    assert result.returncode == 0
    assert os.path.exists(f"{output_dir}/{output_prefix}_interpreted.tsv.tsv")