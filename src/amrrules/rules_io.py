import csv, re
from importlib import resources

def parse_rules_file(rule_file_list):
    # get the correct rules file based on the organism, from the rules directory
    rules_parsed = []
    for rule_file in rule_file_list:
        rule_file_name = f"{rule_file}.tsv"
        try:
            with resources.files("amrrules.rules").joinpath(rule_file_name).open('r', encoding='utf-8') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    rules_parsed.append(row)
        except FileNotFoundError:
            raise FileNotFoundError(f"Rules file '{rule_file_name}' not found in packaged rules/")
    return rules_parsed

def extract_relevant_rules(rules, organism):
    """
    Extract relevant rules for a given organism from the rules list.
    """
    relevant_rules = []
    for rule in rules:
        if rule.get('organism') == organism:
            relevant_rules.append(rule)
    return relevant_rules

def parse_multicopy_rule_mutation(mutation, get_marker=False, gene=None):
    """
    Parse multi-copy rule mutation strings of the form c.[mutation][copy_count], Eg: c.[2611C>T][4]
    Returns a tuple of (base_mutation, copy_count) 
    """

    mutation = mutation.strip()
    mutation_plus_copy = re.match(r"^c\.\[([^\]]+)\]\[(\d+)\]$", mutation)
    if mutation_plus_copy:
        if get_marker:
            marker = gene + f":c.{mutation_plus_copy.group(1)}"
            # just return the marker and the threshold as an int
            return marker, int(mutation_plus_copy.group(2))
        # otherwise return the base mutation as str and the threshold as an int
        return f"c.{mutation_plus_copy.group(1)}", int(mutation_plus_copy.group(2))