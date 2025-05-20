import csv
import os

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

def create_summary_row(sample_rows, sample, rules, no_flag_core):

    summarised_rows = []

    drugs = set() # list of unique drugs to parse
    drug_classes = set() # list of unique drug classes to parse
    for row in sample_rows:
        drug = row.get('drug')
        drug_class = row.get('drug class')
        drugs.add(drug)
        drug_classes.add(drug_class)

    # remove any '-' or '' values
    drugs = [x for x in drugs if x not in ['-', '']]
    drug_classes = [x for x in drug_classes if x not in ['-', '']]
    drugs_and_classes = drugs + drug_classes

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
        matching_rows = [row for row in sample_rows if row.get('drug') == drug_or_class or row.get('drug class') == drug_or_class]
        # if there's only one row, then just add the relevant info to the summary
        if len(matching_rows) == 1:
            # add the sample
            if sample is not None:
                summarised['Name'] = sample
            # get the category, pheno, and evidence grade info
            summarised['category'] = matching_rows[0].get('clinical category')
            summarised['phenotype'] = matching_rows[0].get('phenotype')
            summarised['evidence grade'] = matching_rows[0].get('evidence grade')
            # grab marker and check if it's a core gene, flag if user wants this
            marker_value = matching_rows[0].get('Gene symbol') or matching_rows[0].get('Element symbol')
            if not no_flag_core:
                # check if the gene symbol is a core gene
                context_value = matching_rows[0].get('context')
                if context_value == 'core':
                    marker_value = marker_value + ' (core)'
            summarised['markers'] = marker_value
            summarised['ruleIDs'] = matching_rows[0].get('ruleID')
            summarised['combo rules'] = '-'
            summarised['organism'] = matching_rows[0].get('organism')

            summarised_rows.append(summarised)
        # if there are multiple rows, we need to combine some of the info
        elif len(matching_rows) > 1:
            # add the sample
            if sample is not None:
                summarised['Name'] = sample
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
                if not no_flag_core:
                    if row.get('context') == 'core':
                        gene_symbol = gene_symbol + ' (core)'
                if gene_symbol not in markers:
                    markers.append(gene_symbol)
            summarised['markers'] = ';'.join(markers)
            
            if matched_combo_rules is not None:
                summarised['combo rules'] = ';'.join(combo_ruleIDs)
            else:
                summarised['combo rules'] = '-'
            summarised['organism'] = matching_rows[0].get('organism')
            summarised_rows.append(summarised)
   
    return summarised_rows

def prepare_summary(output_rows, rules, sample_ids, no_flag_core):

    # We now want to write a sumamary file that groups hits based on drug class or drug
    # In each drug/drug class, we want to list the highest category and phenotype for that drug/drug class
    # then all the markers that are associated with that drug/drug class, then all the singleton rule IDs

    # set up the final rows to be output
    summary_rows = []

    # if we've only got one sample, sample_ids will be None
    if sample_ids is None:
        summary_rows = create_summary_row(output_rows, None, rules, no_flag_core)
        return summary_rows

    # we need to process this sample by sample, so the combo rules are correctly assigned
    for sample in sample_ids:
        # extract all the rows for this sample
        sample_rows = [row for row in output_rows if row.get('Name') == sample]
        summarised_row = create_summary_row(sample_rows, sample, rules, no_flag_core)
        summary_rows.append(summarised_row)
    flattened_summary = [item for sublist in summary_rows for item in sublist]
    return flattened_summary

def write_output_files(output_rows, reader, summary_output, args, unmatched_hits, matched_hits):

    # write the output files
    interpreted_output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    #summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')

    minimal_columns = ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade', 'version', 'organism']
    full_columns = ['breakpoint', 'breakpoint standard', 'evidence code', 'evidence limitations', 'PMID', 'rule curation note']
    if args.annot_opts == 'minimal':
        interpreted_output_cols = minimal_columns
    elif args.annot_opts == 'full':
        interpreted_output_cols = minimal_columns + full_columns

    with open(interpreted_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames + interpreted_output_cols, delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} hits did not match a rule.")
    print(f"Output written to {interpreted_output_file}.")

    #with open(summary_output_file, 'w', newline='') as f:
    #    if 'Name' in summary_output[0].keys():
    #        col_names = ['Name', 'drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers', 'ruleIDs', 'combo rules', 'organism']
    #    else:
    #        col_names = ['drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers', 'ruleIDs', 'combo rules', 'organism']
    #    writer = csv.DictWriter(f, fieldnames=col_names, delimiter='\t')
    #    writer.writeheader()
    #    writer.writerows(summary_output)
    #print(f"Summary output written to {summary_output_file}.")