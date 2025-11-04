import re
from amrrules.utils import aa_conversion
from amrrules import __version__

def get_final_rule_matches(matching_rules, type, amrrules_mutation):
    if type == 'Gene presence detected':
        return matching_rules
    elif type == 'Protein variant detected' or type == 'Nucleotide variant detected' or type == "Promoter variant detected":
        final_matching_rules = []
        # now we need to check the mutation, extracting any matching rules
        for rule in matching_rules:
            if rule['mutation'] == amrrules_mutation:
                final_matching_rules.append(rule)
        return final_matching_rules

def check_rules(row, rules, amrfp_nodes):

    # if the row is a point mutation, we need to extract that info and look only for those rules
    # skip any rows that have nothing to do with AMR
    element_type = row.get('Element type') or row.get('Type')
    element_subtype = row.get('Element subtype') or row.get('Subtype')
    if element_type != "AMR":
        return None
    elif element_subtype == "POINT":
        amrrules_mutation, type = extract_mutation(row)
    elif element_subtype == "AMR":
        amrrules_mutation = None
        type = 'Gene presence detected'

    # select rules that match our variation type
    rules_to_check = []
    for rule in rules:
            if rule['variation type'] == type:
                rules_to_check.append(rule)

    # we need to set some logic here about what accessions we're going to be looking for
    # if there's a nodeID, then that's where we want to start
    # otherwise we're going to be checking refseq, or HMM accessions
    #TODO Genbank accessions - not relevant for amrfp as these won't be reported in the output file
    hierarchy_id = row.get('Hierarchy node')
    seq_acc = row.get('Accession of closest sequence') or row.get('Closest reference accession')

    # only grab HMM if it's actually listed
    if row.get('HMM id') != 'NA':
        HMM_acc = row.get('HMM id')
    else:
        HMM_acc = None

    # First we're going to check for the hierarchy ID, and if we have one or matches, we we return that
    matching_rules = [rule for rule in rules_to_check if rule.get('nodeID') == hierarchy_id]
    if len(matching_rules) > 0:
        return(get_final_rule_matches(matching_rules, type, amrrules_mutation))

    # Okay so nothing matched directly to the nodeID, or we would've returned out of the function. 
    # So now we need to check if there's a parent node that matches
    # Note we don't need to check for variation type here, because we're not going to need to go up the hierarchy
    # for that type of rule
    parent_node = amrfp_nodes.get(hierarchy_id)
    while parent_node is not None and parent_node != 'AMR':
        matching_rules = [rule for rule in rules_to_check if rule.get('nodeID') == parent_node]
        if len(matching_rules) > 0:
            return(matching_rules)
        parent_node = amrfp_nodes.get(parent_node)

    #Okay so using the nodeID didn't work, so now we need to check the sequence accession
    # start with the nucleotide accessions
    matching_rules = [rule for rule in rules_to_check if rule.get('nucleotide accession') == seq_acc]
    if len(matching_rules) > 0:
        return(get_final_rule_matches(matching_rules, type, amrrules_mutation))
    # then check the protein accessions
    matching_rules = [rule for rule in rules_to_check if rule.get('protein accession') == seq_acc]
    if len(matching_rules) > 0:
        return(get_final_rule_matches(matching_rules, type, amrrules_mutation))

    #TODO: HMM accession check

    # if nothing matched, then we reurn None
    return None

def annotate_rule(row, rules, annot_opts, version=__version__):
    minimal_columns = ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade', 'version', 'organism']
    full_columns = ['breakpoint', 'breakpoint standard', 'breakpoint condition', 'evidence code', 'evidence limitations', 'PMID', 'rule curation note']

    if rules is None:
        # if we didn't find a matching rule, then we need to add new columns for each of the options but using '-' as the value
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = '-'
        elif annot_opts == 'full':
            for col in minimal_columns + full_columns:
                row[col] = '-'
        row['version'] = version
        row['organism'] = '-'
        return [row]
    # if we found multiple rules, we want to return each rule as its own row in the output file
    if len(rules) > 1:
        # we need to create a new row for each rule
        output_rows = []
        for rule in rules:
            new_row = row.copy()
            if annot_opts == 'minimal':
                for col in minimal_columns:
                    new_row[col] = rule.get(col)
            elif annot_opts == 'full':
                for col in minimal_columns + full_columns:
                    new_row[col] = rule.get(col)
            new_row['version'] = version
            output_rows.append(new_row)
        return output_rows
    # if we found a single rule, we want to annotate the row with the rule info
    else:
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = rules[0].get(col)
        elif annot_opts == 'full':
            for col in minimal_columns + full_columns:
                row[col] = rules[0].get(col)
        row['version'] = version
        return [row]