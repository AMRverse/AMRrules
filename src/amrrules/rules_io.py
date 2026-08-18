from collections import defaultdict
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

def extract_unknown_core_rules(rules, card_drug_map):
    """
    Extract rules that are for core resistance mechanisms that have unknown genes.
    Group the rules by drug class, then drug, as we have done within each sample.
    """
    unknown_rules = defaultdict(lambda: defaultdict(list))
    for rule in rules:
        if rule.get('gene') in ('unknown', 'none') and rule.get('phenotype') == 'wildtype':
            if rule.get('drug class') == '-':
                rule['drug class'] = card_drug_map.get(rule.get('drug'))
            unknown_rules[rule.get('drug class')][rule.get('drug')].append(rule)
    return unknown_rules

def parse_multicopy_rule_mutation(mutation, get_marker=False, gene=None):
    """
    Parse multi-copy rule mutation strings. Two formats are supported:
      - c.[mutation][copy_count]  e.g. c.[2611C>T][4]
        -> a specific point mutation, required for "Nucleotide variant detected in multi-copy gene"
      - c.[copy_count]            e.g. c.[2]
        -> the gene itself present in >= copy_count copies, independent of
           any specific mutation ("Gene copy number variant detected")
    Returns a tuple of (base_mutation_or_marker, copy_count), or
    (None, None) if the string doesn't match either format.
    """

    mutation = mutation.strip()

    # mutation variant form: c.[mutation][copy_count]
    mutation_plus_copy = re.match(r"^c\.\[([^\]]+)\]\[(\d+)\]$", mutation)
    if mutation_plus_copy:
        if get_marker:
            marker = gene + f":c.{mutation_plus_copy.group(1)}"
            # just return the marker and the threshold as an int
            return marker, int(mutation_plus_copy.group(2))
        # otherwise return the base mutation as str and the threshold as an int
        return f"c.{mutation_plus_copy.group(1)}", int(mutation_plus_copy.group(2))

    # gene presence form: c.[copy_count]
    copy_only = re.match(r"^c\.\[(\d+)\]$", mutation)
    if copy_only:
        if get_marker:
            # the marker for a whole-gene copy-number rule is just the gene itself
            return gene, int(copy_only.group(1))
        return None, int(copy_only.group(1))

    # not a multicopy-format mutation string at all
    return None, None

def get_combination_rules(rules, organism):
    """
    Extracts combination rules for a given organism.
    """
    combo_rules = [
        r for r in rules
        if r.get('organism') == organism
        and r.get('variation type') in ('Combination')
    ]
    return combo_rules

def evaluate_logic_string(logic_string, id_list):
    """
    Evaluates a logic string against a list of IDs using strict word boundaries.

    Args:
        logic_string (str): A string containing logical expressions (e.g., "ECO1016 & ECO1026").
        id_list (list/set): A collection of IDs to compare against.

    Returns:
        bool: True if the logic evaluates to True, False otherwise.
    """
    id_set = set(id_list)

    # 1. Convert logical operators to Python equivalents
    # Use word boundaries for 'AND'/'OR' to avoid messing up IDs containing 'AND' or 'OR'
    python_logic = logic_string.replace('&', ' and ').replace('|', ' or ')

    # 2. Extract all distinct alphanumeric tokens (IDs) from the logic string
    # This automatically ignores parentheses, spaces, and operators
    tokens_in_logic = set(re.findall(r'\b\w+\b', logic_string))

    # 3. Safely substitute each ID with its membership check
    for id_ in tokens_in_logic:
        # Skip Python keywords generated from operators
        if id_ in ('and', 'or', 'not'):
            continue
            
        # \b ensures exact match (e.g., matches "NGO006" but NOT "NGO0065")
        pattern = r'\b' + re.escape(id_) + r'\b'
        replacement = f"('{id_}' in id_set)"
        python_logic = re.sub(pattern, replacement, python_logic)

    # 4. Safely evaluate the expression
    try:
        return eval(python_logic, {"__builtins__": None}, {"id_set": id_set})
    except Exception as e:
        raise ValueError(f"Error evaluating logic string: {logic_string}") from e