import os
import shutil
import csv

SRC_DIR = "rules"
DST_DIR = "src/amrrules/rules"

def copy_rules():
    if not os.path.isdir(SRC_DIR):
        raise FileNotFoundError(f"Source rules folder not found: {SRC_DIR}")

    os.makedirs(DST_DIR, exist_ok=True)

    for file in os.listdir(SRC_DIR):
        if file.endswith(".txt"):
            shutil.copy(os.path.join(SRC_DIR, file), os.path.join(DST_DIR, file))
            print(f"Copied: {file}")

def create_rules_file_key():
    """
    Create a key file that lists, for each unique organism, the rules file name.
    The file will look like this:
    s__Klebsiella pneumoniae    Klebsiella_pneumoniae.txt
    s__Escherichia coli        Escherichia_coli.txt
    """

    if DST_DIR is None:
        # Default path is relative to this file (src/AMRrules/rules/)
        rule_dir = os.path.join(os.path.dirname(__file__), "rules")
    
    if not os.path.isdir(DST_DIR):
        raise FileNotFoundError(f"Rules directory not found: {DST_DIR}")
    
    key_file_rows = {}

    for filename in os.listdir(DST_DIR):
        if filename.endswith(".txt"):
            file_org_name = os.path.splitext(filename)[0]
            #extract all unique values from the 'orgnanism' column
            reader = csv.DictReader(open(os.path.join(DST_DIR, filename), 'r'), delimiter='\t')
            organisms = set()
            for row in reader:
                organisms.add(row.get('organism'))
            # now add each organism as a row to the key file
            for organism in organisms:
                    key_file_rows[organism] = file_org_name
    # write the key file
    key_file_path = os.path.join(DST_DIR, "rule_key_file.tsv")
    with open(key_file_path, 'w') as key_file:
        writer = csv.writer(key_file, delimiter='\t')
        for organism, rules_file in key_file_rows.items():
            writer.writerow([organism, rules_file])

    print(f"Rules key file created at: {key_file_path}")

if __name__ == "__main__":
    copy_rules()
    create_rules_file_key()
