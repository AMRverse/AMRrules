def check_comboo_rules(ruleIDs, rules, drug, drug_class):

    matched_combo_rules = []
    # extract all the combination rules from the rules file for this drug or drug class
    if drug is not None:
        combo_rules = [rule for rule in rules if rule.get('drug') == drug and rule.get('variation type') == 'Combination']
    elif drug_class is not None:
        combo_rules = [rule for rule in rules if rule.get('drug class') == drug_class and rule.get('variation type') == 'Combination']
    
    # for each rule, we need to extract the ruleID logic and check if it matches our ruleIDs
    for rule in combo_rules:
        ruleID_logic = rule.get('gene')
        # ruleID logic is a string of the form "gene1 & gene2 | gene3"
        # we need to replace the & with a python 'and' and the | with a python 'or'
        matched_combo = evaluate_logic_string(ruleID_logic, ruleIDs)
        if matched_combo:
            # add to the list of matched combos
            matched_combo_rules.append(ruleID_logic)
    # if we didn't find any matching combo rules, then we return None
    if len(matched_combo_rules) == 0:
        return None
    else:
        return matched_combo_rules

def evaluate_logic_string(logic_string, id_list):
    """
    Evaluates a logic string against a list of IDs.

    Args:
        logic_string (str): A string containing logical expressions (e.g., "ECO1016 & ECO1026").
        id_list (list): A list of IDs to compare against (e.g., ["ECO1016", "ECO1026"]).

    Returns:
        bool: True if the logic evaluates to True, False otherwise.
    """
    # Convert the list of IDs to a set for efficient membership testing
    id_set = set(id_list)

    # Replace logical operators in the string with Python equivalents
    python_logic = logic_string.replace('&', ' and ').replace('|', ' or ')

    # Wrap each ID in the logic string with a check for membership in the ID set
    for id_ in id_set.union(set(logic_string.replace('(', '').replace(')', '').split())):
        python_logic = python_logic.replace(id_, f"'{id_}' in id_set")

    # Evaluate the logic string
    try:
        return eval(python_logic)
    except Exception as e:
        raise ValueError(f"Error evaluating logic string: {logic_string}") from e

def write_summary(output_rows, rules):

    # We now want to write a sumamary file that groups hits based on drug class or drug
    # In each drug/drug class, we want to list the highest category and phenotype for that drug/drug class
    # then all the markers that are associated with that drug/drug class, then all the singleton rule IDs
    #TODO: implement combo rules

    drugs = [] # list of unique drugs to parse
    drug_classes = [] # list of unique drug classes to parse
    for row in output_rows:
        drug = row.get('drug')
        drug_class = row.get('drug class')
        drugs.append(drug)
        drug_classes.append(drug_class)
    # remove duplicates
    drugs = list(set(drugs))
    drug_classes = list(set(drug_classes))
    # remove any '-' or '' values
    drugs = [x for x in drugs if x not in ['-', '']]
    drug_classes = [x for x in drug_classes if x not in ['-', '']]
    drugs_and_classes = drugs + drug_classes

    summary_rows = []
    for drug_or_class in drugs_and_classes:
        if drug_or_class in drugs:
            summarised = {'drug': drug_or_class, 'drug class': '-'} #TODO: look up drug class for drug from the card ontology
            drug = drug_or_class
            drug_class = None
        elif drug_or_class in drug_classes:
            summarised = {'drug': '-', 'drug class': drug_or_class}
            drug = None
            drug_class = drug_or_class
        # extract all rows that match this drug or class
        matching_rows = [row for row in output_rows if row.get('drug') == drug_or_class or row.get('drug class') == drug_or_class]
        # if there's only one row, then just add the relevant info to the summary
        if len(matching_rows) == 1:
            summarised['category'] = matching_rows[0].get('clinical category')
            summarised['phenotype'] = matching_rows[0].get('phenotype')
            summarised['evidence grade'] = matching_rows[0].get('evidence grade')
            summarised['markers'] = matching_rows[0].get('Gene symbol') or matching_rows[0].get('Element symbol')
            summarised['ruleIDs'] = matching_rows[0].get('ruleID')
            summary_rows.append(summarised)
        # if there are multiple rows, we need to combine some of the info
        elif len(matching_rows) > 1:
            # combine all the rule IDs into a single string
            ruleIDs = []
            for row in matching_rows:
                ruleID = row.get('ruleID')
                if ruleID not in ruleIDs:
                    ruleIDs.append(ruleID)
            summarised['ruleIDs'] = ';'.join(ruleIDs)

            # we now need to check if there are any combo rules, and if so, we need to add them to the summary
            matched_combo_rules = check_comboo_rules(ruleIDs, rules, drug, drug_class)

            # extract the current categories, phenotypes and evidence grades
            categories = [row.get('clinical category') for row in matching_rows]
            phenotypes = [row.get('phenotype') for row in matching_rows]
            evidence_grades = [row.get('evidence grade') for row in matching_rows]

            # if there are combo rules, then we need to add these rules to our clinical category and phenotype checks
            # and return the highest category and phenotype including the combos
            if matched_combo_rules is not None:
                combo_ruleIDs = []
                # we need to check if the combo rules are in the matching rows
                for rule in matched_combo_rules:
                    # extract the category, phenotype and evidence grade from the rule
                    for r in rules:
                        if r.get('gene') == rule:
                            categories.append(r.get('clinical category'))
                            phenotypes.append(r.get('phenotype'))
                            evidence_grades.append(r.get('evidence grade'))
                            # extract the ruleID that links to this combo rule (rather than specifying the actual combo rule)
                            combo_ruleIDs.append(r.get('ruleID'))
    
            # get the highest category
            highest_category = max(categories, key=lambda x: ['-', 'S', 'I', 'R'].index(x))
            summarised['category'] = highest_category
            # get the highest phenotype
            highest_phenotype = max(phenotypes, key=lambda x: ['-', 'wildtype', 'nonwildtype'].index(x))
            summarised['phenotype'] = highest_phenotype
            # get the highest evidence grade
            highest_evidence_grade = max(evidence_grades, key=lambda x: ['-', 'weak', 'moderate', 'strong'].index(x))
            summarised['evidence grade'] = highest_evidence_grade
            # combine all the markers into a single string
            markers = []
            for row in matching_rows:
                gene_symbol = row.get('Gene symbol') or row.get('Element symbol')
                if gene_symbol not in markers:
                    markers.append(gene_symbol)
            summarised['markers'] = ';'.join(markers)
            
            if matched_combo_rules is not None:
                summarised['combo rules'] = ';'.join(combo_ruleIDs)
            else:
                summarised['combo rules'] = '-'
            summary_rows.append(summarised)
    
    return summary_rows