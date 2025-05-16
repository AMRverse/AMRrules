import csv
import urllib.request
from importlib import resources
import io

def download_and_parse_reference_gene_hierarchy(url):
    print(f"Downloading Reference Gene Hierarchy from {url}...")
    
    # Use urllib to download the file
    with urllib.request.urlopen(url) as response:
        content = response.read().decode('utf-8')  # Decode the response as UTF-8
    
    # Parse the downloaded content as a TSV file
    content_io = io.StringIO(content)
    reader = csv.DictReader(content_io, delimiter='\t')
    
    amrfp_nodes = {}
    for row in reader:
        node_id = row.get('node_id')
        parent_node = row.get('parent_node_id')
        amrfp_nodes[node_id] = parent_node

    print(f"Downloaded and parsed {len(amrfp_nodes)} rows from the Reference Gene Hierarchy.")
    return amrfp_nodes

def parse_rules_file(rule_file_list):
    # get the correct rules file based on the organism, from the rules directory
    rules_parsed = []
    for rule_file in rule_file_list:
        rule_file_name = f"{rule_file}.txt"
        try:
            with resources.files("amrrules.rules").joinpath(rule_file_name).open('r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    rules_parsed.append(row)
        except FileNotFoundError:
            raise FileNotFoundError(f"Rules file '{rule_file_name}' not found in packaged rules/")
    return rules_parsed
    
def get_organisms(organism_file):

    organism_dict = {}
    with open(organism_file, 'r') as org_file:
        for line in org_file:
            sample_id, organism_name = line.strip().split('\t')
            if sample_id not in organism_dict:
                organism_dict[sample_id] = organism_name
            elif sample_id in organism_dict:
                raise ValueError(f"Duplicate sample ID found in organism file: {sample_id}. Please ensure that each sample ID is unique.")
    return organism_dict

def extract_relevant_rules(rules, organism):
    """
    Extract relevant rules for a given organism from the rules list.
    """
    relevant_rules = []
    for rule in rules:
        if rule.get('organism') == organism:
            relevant_rules.append(rule)
    return relevant_rules