# This script extracts the names of organisms from text files in the rules directory
# and writes them to a Python file (utils.py) in the src/amrrules directory.
# this script is run during the build process to ensure that the utils.py file is up to date

import os

def extract_organisms(rules_dir, output_file):
    # Get all .txt files in the rules directory
    filenames = [
        os.path.splitext(f)[0]
        for f in os.listdir(rules_dir)
        if f.endswith(".txt")
    ]

    # Write the filenames to utils.py
    with open(output_file, "w") as f:
        f.write("# This file is auto-generated during the build process\n")
        f.write("allowed_organisms = ")
        f.write(repr(filenames))
        f.write("\n")

if __name__ == "__main__":
    # Define the rules directory and output file
    rules_dir = "rules"  # Adjust this path if necessary
    output_file = "src/amrrules/utils.py"
    extract_organisms(rules_dir, output_file)