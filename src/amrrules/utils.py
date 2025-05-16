import os, csv
from importlib import resources
import warnings

aa_conversion = {'G': 'Gly', 'A': 'Ala', 'S': 'Ser', 'P': 'Pro', 'T': 'Thr', 'C': 'Cys', 'V': 'Val', 'L': 'Leu', 'I': 'Ile', 
                 'M': 'Met', 'N': 'Asn', 'Q': 'Gln', 'K': 'Lys', 'R': 'Arg', 'H': 'His', 'D': 'Asp', 'E': 'Glu', 'W': 'Trp', 
                 'Y': 'Tyr', 'F': 'Phe', '*': 'STOP'}

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
    Read the organism file and return a dictionary of sample IDs and their corresponding organism names."""

    organism_dict = {}
    with open(organism_file, 'r') as org_file:
        for line in org_file:
            sample_id, organism_name = line.strip().split('\t')
            if sample_id not in organism_dict:
                organism_dict[sample_id] = organism_name
            elif sample_id in organism_dict:
                raise ValueError(f"Duplicate sample ID found in organism file: {sample_id}. Please ensure that each sample ID is unique.")
    return organism_dict

def validate_input_file(input_file, organism_file, organism_user, tool='amrfp'):
    """
    Validate the input file to ensure it contains the required columns.
    """
    #TODO: these checks will need to change if the tool isn't amrfp
    with open(input_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError("Input file does not contain 'Hierarchy node' column. Please re-run AMRFinderPlus with the --print_node option to ensure this column is in the output file.")
        if organism_file and 'Name' not in reader.fieldnames:
            raise ValueError("Input file does not contain 'Name' column, which is required when using an organism file. Please ensure the input file has a 'Name' column, and that the IDs in this column match the IDs in column 1 of your organism file.")

        # extract the sample IDs from the intput file
        if 'Name' in reader.fieldnames:
            sample_ids = set()
            for row in reader:
                sample_ids.add(row.get('Name'))

        # if we have an organism file, we need to create a dictionary (which we will then pass out to the main engine), and check all those IDs are in our organism file
        if organism_file:
            organism_dict = get_organisms(organism_file)
            # check that all the sample IDs have a corresponding organism
            # for now raise an error if any are missing
            samples_to_skip = set()
            for sample_id in sample_ids:
                if sample_id not in organism_dict.keys():
                    samples_to_skip.add(sample_id)
            if len(samples_to_skip) > 0:
                warnings.warn(f"Sample IDs '{', '.join(samples_to_skip)}' from input file not found in organism file. These samples will be skipped for interpretation.")
                # create a new list of samples that will be used for interpretation
                sample_ids = sample_ids - samples_to_skip
            else:
                print(f"All sample IDs from input file found in organism file. Proceeding with interpretation.")
        # if the user has supplied just organism, then we don't need to check for samples, and we can just create a dummy dictionary
        # for picking what rules to use
        if organism_user:
            organism_dict = {'': organism_user}
            sample_ids = None
       
    return organism_dict, sample_ids