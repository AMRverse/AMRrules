import os

aa_conversion = {'G': 'Gly', 'A': 'Ala', 'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'C': 'Cys', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 
                 'M': 'Met', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys', 'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu', 'W': 'Trp', 
                 'Y': 'Tyr', 'F': 'Phe', '*': 'STOP'}

def get_supported_organisms(rule_dir: str = None):
    """
    Return a list of organism names by scanning the rules folder.
    Assumes filenames are of the form 'organism.txt'.
    """
    if rule_dir is None:
        # Default path is relative to this file (src/AMRrules/rules/)
        rule_dir = os.path.join(os.path.dirname(__file__), "rules")
    
    if not os.path.isdir(rule_dir):
        raise FileNotFoundError(f"Rules directory not found: {rule_dir}")

    organisms = []
    for filename in os.listdir(rule_dir):
        if filename.endswith(".txt"):
            organisms.append(os.path.splitext(filename)[0])  # Strip .txt
    return sorted(organisms)