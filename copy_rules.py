import os
import shutil

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

if __name__ == "__main__":
    copy_rules()
