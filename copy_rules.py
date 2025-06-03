import os
import shutil
import csv

SRC_DIR = "rules"
DST_DIR = "src/amrrules/rules"

def clean_and_copy_rules():
    # run some checks on each rules file before it's copied
    # make sure that there are no trailing spaces in any of the columns
    # make sure there are no empty rows
    # make sure that the second row isn't the 'required' 'optional' stuff from the spec
    if not os.path.isdir(SRC_DIR):
        raise FileNotFoundError(f"Source rules folder not found: {SRC_DIR}")

    os.makedirs(DST_DIR, exist_ok=True)

    for file in os.listdir(SRC_DIR):
        if file.endswith(".txt"):
            file_path = os.path.join(SRC_DIR, file)
            out_file_path = os.path.join(DST_DIR, file)
            with open(file_path, 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                cleaned_rows = []
                for row in reader:
                    # skip empty rows
                    if not row or all(cell.strip() == '' for cell in row):
                        continue
                    # skip the second row if it contains 'required' or 'optional'
                    if len(row) > 1 and (row[0].lower() == 'required' or row[0].lower() == 'optional'):
                        continue
                    # strip trailing spaces from each cell
                    cleaned_row = [cell.strip() for cell in row]
                    cleaned_rows.append(cleaned_row)
            # write the cleaned rows back to the file
            with open(out_file_path, 'w', newline='') as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerows(cleaned_rows)
            print(f"Cleaned and copied: {file}")

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
    clean_and_copy_rules()
    create_rules_file_key()
