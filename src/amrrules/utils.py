import os, csv
from importlib import resources
import warnings

aa_conversion = {'G': 'Gly', 'A': 'Ala', 'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'C': 'Cys', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 
                 'M': 'Met', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys', 'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu', 'W': 'Trp', 
                 'Y': 'Tyr', 'F': 'Phe', '*': 'STOP'}

minimal_columns = ['ruleID', 'gene context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade', 'version', 'organism']
full_columns = ['breakpoint', 'breakpoint standard', 'breakpoint condition', 'evidence code', 'evidence limitations', 'PMID', 'rule curation note']

CATEGORY_ORDER = ['S', 'I', 'R']
PHENOTYPE_ORDER = ['wildtype', 'nonwildtype']
EVIDENCE_GRADE_ORDER = ['very low', 'low', 'moderate', 'high']

def get_supported_organisms(rule_dir: str = None):
    """
    Return a list of organism names by scanning organism names in the rules folder.
    """
    rule_dir = resources.files("amrrules.rules")
    
    if rule_dir is None:
        raise ValueError("Rules directory not found. Please check your installation.")

    organisms = set()
    for entry in rule_dir.iterdir():
        if entry.name != "rule_key_file.tsv" and entry.name.endswith(".txt"):
            # Extract organism name from filename
            reader = csv.DictReader(open(entry, 'r'), delimiter='\t')
            for row in reader:
                organisms.add(row.get('organism'))
    return sorted(organisms)

def get_organisms(organism_file):
    """
    Read the organism file and return a dictionary of sample IDs and their corresponding organism names.
    """

    organism_dict = {}
    with open(organism_file, 'r') as org_file:
        for line in org_file:
            sample_id, organism_name = line.strip().split('\t')
            if sample_id not in organism_dict:
                organism_dict[sample_id] = organism_name
            elif sample_id in organism_dict:
                raise ValueError(f"Duplicate sample ID found in organism file: {sample_id}. Please ensure that each sample ID is unique.")
    return organism_dict

def validate_amrfp_file(amrfp_file, multi_entry=False):
    """
    Validate that the AMRFinderPlus input file contains the required columns. All files must have Hierarchy node. Multi entry files must have Name column.
    """
    with open(amrfp_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        samples_in_file = '' # set this to empty string if we only have one sample and no 'Name' column
        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError(f"Input AMRFinderPlus file is missing required column: 'Hierarchy node'. Please re-run AMRFinderPlus with the --print_node option to ensure this column is in the output file.")
        if multi_entry and 'Name' not in reader.fieldnames:
            raise ValueError(f"Input AMRFinderPlus file is missing required column: 'Name'. Please ensure this column is present so we can match samples to organisms in the supplied organism file.")
        if 'Name' in reader.fieldnames:
            samples_in_file = set()
            for row in reader:
                samples_in_file.add(row.get('Name'))
    return samples_in_file

def check_sample_ids(samples_with_org, samples_in_input):
    """
    Check that all samples in the input file have a corresponding organism in the organism file.
    It's okay if there are samples in the organism file that aren't in the input file, we can just raise a warning in that instance.
    """
    
    missing_samples = samples_in_input - samples_with_org
    if missing_samples:
        raise ValueError(f"The following sample IDs from the input file are missing in the organism file: {', '.join(missing_samples)}. Please ensure all sample IDs are present in the organism file.")
    missing_from_input = samples_with_org - samples_in_input
    if missing_from_input:
        warnings.warn(f"The following sample IDs from the organism file are not present in the input file:\n{'\n'.join(missing_from_input)}\nAs there are no entries in the input file for these samples, they won't have interpretation results. Please check your input file if this is not what you expect.")
    return True