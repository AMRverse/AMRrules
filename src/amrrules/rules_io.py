import csv
from importlib import resources

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

def extract_relevant_rules(rules, organism):
    """
    Extract relevant rules for a given organism from the rules list.
    """
    relevant_rules = []
    for rule in rules:
        if rule.get('organism') == organism:
            relevant_rules.append(rule)
    return relevant_rules