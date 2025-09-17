import csv
import os
from amrrules.resources import ResourceManager as rm

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

def summarise_by_drug_or_class(row_indicies, row, drug_class_name):
    
    if drug_class_name not in row_indicies:
            row_indicies[drug_class_name] = [row]
    else:
        row_indicies[drug_class_name].append(row)
    return row_indicies

def create_summary_row(sample_rows, sample, rules, no_flag_core, no_rule_interpretation):

    summarised_rows = []

    drugs = set() # list of unique drugs to parse
    drug_classes = set() # list of unique drug classes to parse

    card_amrfp_conversion = rm().get_amrfp_card_conversion() # get the AMRFP to CARD conversion mapping
    card_drug_map = rm().get_card_drug_class_map() # get CARD drugs and their associated classes

    # we need to evaluate markers by drug/class, because that's how we're summarising
    # so first every marker/row needs to be assigned to a drug or class
    # for rows where rules have been applied, that's easy, it's already done
    # for rows where they haven't, we need to convert the AMRFP Subclass to it's corresponding CARD drug or class
    # some rows will have NA for Subclass because they aren't relevant, these can be skipped
    # some rows will have a Subclass but be NA for card, these need to be grouped under 'other markers'

    drug_class_row_indicies = {}

    for row in sample_rows:
        if row.get('ruleID') == '-' and row.get('Subclass') != 'NA': # convert AMRFP Subclass to CARD, skipping non-AMR rows
            subclasses = row.get('Subclass')
            drug_class = None
            # if there are '/' between items, we need to split these
            subclasses = subclasses.split('/') if '/' in subclasses else [subclasses]
            for subclass in subclasses:
                summary_drug = card_amrfp_conversion.get(subclass, {}).get('drug', 'NA')
                if summary_drug != 'NA' and summary_drug != '-':
                    drugs.add(summary_drug)
                if summary_drug == 'NA' or summary_drug == '-': # only grab class if we can't get a drug
                    drug_class = card_amrfp_conversion.get(subclass, {}).get('class', 'NA')
                    if summary_drug == 'NA' and drug_class == 'NA': # both are NA so it's other
                        summary_drug = 'other markers'
                    elif summary_drug != 'NA' or summary_drug != '-': # only drug is NA, so we can set the class instead
                        summary_drug = drug_class
                        drug_classes.add(drug_class)
                # because no rule has been applied, we need to fill in some of the additional info
                # this should be set by the user
                if no_rule_interpretation == 'nwtR':
                    row['phenotype'] = 'nonwildtype'
                    row['clinical category'] = 'R'
                    row['evidence grade'] = 'very low'
                elif no_rule_interpretation == 'nwtS':
                    row['phenotype'] = 'nonwildtype'
                    row['clinical category'] = 'S'
                    row['evidence grade'] = 'very low'
                drug_class_row_indicies = summarise_by_drug_or_class(drug_class_row_indicies, row, summary_drug)
        elif row.get('ruleID') != '-': # process the rows with matching rules
            summary_drug = row.get('drug')
            if summary_drug == '-':
                summary_drug = row.get('drug class')
                drug_classes.add(summary_drug)
            else:
                drugs.add(summary_drug)
            drug_class_row_indicies = summarise_by_drug_or_class(drug_class_row_indicies, row, summary_drug)
        else: # these are rows we want to skip
            continue

    # once each row is assigned to a drug or class, we can then process all rows for that drug/class
    for drug_or_class in drug_class_row_indicies.keys():
        # initialise values
        all_markers_no_rule = False
        drug = None
        drug_class = None
        # work out if we are dealing with a drug or a class
        if drug_or_class in drugs:
            drug = drug_or_class
            drug_class = card_drug_map.get(drug_or_class, '-')
            summarised = {'drug': drug_or_class, 'drug class': drug_class}
        elif drug_or_class in drug_classes:
            summarised = {'drug': '-', 'drug class': drug_or_class}
            drug = None
            drug_class = drug_or_class
        elif drug_or_class == 'other markers':
            summarised = {'drug': '-', 'drug class': 'other markers'}
            # no need to apply combo rules in this case
            drug = None
            drug_class = None 
        rows_to_process = drug_class_row_indicies[drug_or_class]
        # if there's only one row, then just add the relevant info to the summary
        if len(rows_to_process) == 1:
            # add the sample
            if sample is not None:
                summarised['Name'] = sample
            # get the category, pheno, and evidence grade info
            summarised['category'] = rows_to_process[0].get('clinical category')
            phenotype = rows_to_process[0].get('phenotype')
            summarised['phenotype'] = phenotype
            summarised['evidence grade'] = rows_to_process[0].get('evidence grade')
            # grab marker and check if it's a core gene, flag if user wants this
            marker_value = rows_to_process[0].get('Gene symbol') or rows_to_process[0].get('Element symbol')
            context_value = rows_to_process[0].get('context')
            if context_value == 'core' and not no_flag_core:
                marker_value = marker_value + ' (core)'
            # if the marker is wildtype, add it to the wt marker column
            if phenotype == 'wildtype':
                summarised['wt markers'] = marker_value
                summarised['markers (with rule)'] = '-'
                summarised['markers (no rule)'] = '-'
            # if it's not, put the marker in either a rule or no rule column
            # depending on if there's a rule
            ruleID = rows_to_process[0].get('ruleID')
            if phenotype == 'nonwildtype':
                summarised['wt markers'] = '-'
                if ruleID == '-':
                    summarised['markers (no rule)'] = marker_value
                    summarised['markers (with rule)'] = '-'
                else:
                    summarised['markers (with rule)'] = marker_value
                    summarised['markers (no rule)'] = '-'

            summarised['ruleIDs'] = ruleID
            summarised['combo rules'] = '-'
            summarised['organism'] = rows_to_process[0].get('organism')

            summarised_rows.append(summarised)
        # if there are multiple rows, we need to combine some of the info
        # note that there may be some rows where rules have been applied, and some where they haven't
        elif len(rows_to_process) > 1:
            # add the sample
            if sample is not None:
                summarised['Name'] = sample
            # combine all the rule IDs into a single string
            ruleIDs = []
            for row in rows_to_process:
                ruleID = row.get('ruleID')
                # then all markers have no rules
                if ruleID == '-':
                    all_markers_no_rule = True
                if ruleID not in ruleIDs and ruleID != '-':
                    ruleIDs.append(ruleID)
            if all_markers_no_rule: # want to make sure this is at the end of the string
                summarised['ruleIDs'] = '-'
            elif len(ruleIDs) == 1:
                summarised['ruleIDs'] = ruleIDs[0]
            else:
                summarised['ruleIDs'] = ';'.join(ruleIDs)

            # we now need to check if there are any combo rules, and if so, we need to add them to the summary
            if drug or drug_class:
                matched_combo_rules = check_comboo_rules(ruleIDs, rules, drug, drug_class)
            else:
                matched_combo_rules = None

            # extract the current categories, phenotypes and evidence grades
            categories = [row.get('clinical category') for row in rows_to_process]
            phenotypes = [row.get('phenotype') for row in rows_to_process]
            evidence_grades = [row.get('evidence grade') for row in rows_to_process]

            # if there are combo rules, then we need to add these rules to our clinical category and phenotype checks
            # and return the highest category and phenotype including the combos
            if matched_combo_rules:
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
            highest_evidence_grade = max(evidence_grades, key=lambda x: ['-', 'very low', 'low', 'moderate', 'high'].index(x))
            summarised['evidence grade'] = highest_evidence_grade
            
            # combine all the markers into three strings
            # one for wt markers
            # one for markers with a rule
            # one for markers without a rule
            wt_markers = []
            markers_rule = []
            markers_no_rule = []
            for row in rows_to_process:
                gene_symbol = row.get('Gene symbol') or row.get('Element symbol')
                if not no_flag_core:
                    if row.get('context') == 'core':
                        gene_symbol = gene_symbol + ' (core)'
                if row.get('phenotype') == 'wildtype':
                    if gene_symbol not in wt_markers:
                        wt_markers.append(gene_symbol)
                elif row.get('phenotype') == 'nonwildtype':
                    if row.get('ruleID') == '-':
                        if gene_symbol not in markers_no_rule:
                            markers_no_rule.append(gene_symbol)
                    else:
                        if gene_symbol not in markers_rule:
                            markers_rule.append(gene_symbol)
            # set the lists to blank if there are no entries
            if wt_markers == []:
                wt_markers = ['-']
            if markers_rule == []:
                markers_rule = ['-']
            if markers_no_rule == []:
                markers_no_rule = ['-']
            summarised['wt markers'] = ';'.join(wt_markers)
            summarised['markers (with rule)'] = ';'.join(markers_rule)
            summarised['markers (no rule)'] = ';'.join(markers_no_rule)

            if matched_combo_rules is not None:
                summarised['combo rules'] = ';'.join(combo_ruleIDs)
            else:
                summarised['combo rules'] = '-'
            summarised['organism'] = rows_to_process[0].get('organism')
            summarised_rows.append(summarised)
   
    return summarised_rows

def order_rows(rows):

    # First, we want to sort by drug class, ensuring that 'other markers' is at the end
    rows = sorted(rows, key=lambda x: x.get('drug class') == 'other markers')

    # Then, we want to group by drug class
    grouped = {}
    for row in rows:
        drug_class = row.get('drug class')
        if drug_class not in grouped:
            grouped[drug_class] = []
        grouped[drug_class].append(row)

    # Now we want to order each group by drug
    # when we group by drug, we want to put the '-' option at the end
    for drug_class in grouped:
        grouped[drug_class] = sorted(grouped[drug_class], key=lambda x: x.get('drug') == '-')

    # now we need to flatten this back out, so that we remove the upper level key
    flattened = [item for sublist in grouped.values() for item in sublist]

    return flattened

def prepare_summary(output_rows, rules, sample_ids, no_flag_core, no_rule_interpretation):

    # We now want to write a sumamary file that groups hits based on drug class or drug
    # In each drug/drug class, we want to list the highest category and phenotype for that drug/drug class
    # then all the markers that are associated with that drug/drug class, then all the singleton rule IDs

    # set up the final rows to be output
    summary_rows = []

    # if we've only got one sample, sample_ids will be None
    if sample_ids is None:
        summary_rows = create_summary_row(output_rows, None, rules, no_flag_core, no_rule_interpretation)
        ordered_summary_rows = order_rows(summary_rows)
        #flattened_summary = [item for sublist in ordered_summary_rows for item in sublist]
        return ordered_summary_rows

    # we need to process this sample by sample, so the combo rules are correctly assigned
    for sample in sample_ids:
        # extract all the rows for this sample
        sample_rows = [row for row in output_rows if row.get('Name') == sample]
        summarised_row = create_summary_row(sample_rows, sample, rules, no_flag_core)
        # order the rows for this sample by drug and class
        ordered_summarised_row = order_rows(summarised_row)
        summary_rows.append(ordered_summarised_row)
    flattened_summary = [item for sublist in summary_rows for item in sublist]
    return flattened_summary