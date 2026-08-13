from collections import defaultdict
import re
from amrrules.rules_io import parse_multicopy_rule_mutation, get_combination_rules, evaluate_logic_string
from amrrules.genotype_parser import Genotype


def _best_tier_rule(candidates, observed_copies):
    """
    candidates: list of (threshold, rule) tuples for the same family.
    Returns the rule with the highest threshold <= observed_copies, or
    None if no rule qualifies (e.g. only 1 copy observed, but the lowest
    defined tier requires 2).
    """
    qualifying = [(t, r) for t, r in candidates if t is not None and t <= observed_copies]
    if not qualifying:
        return None
    return max(qualifying, key=lambda tr: tr[0])[1]


def apply_copy_number_rules(geno_objs, rules, card_drug_map):
    """
    geno_objs: list of Genotype objects for a sample
    rules: the complete rules list for the relevant organism, not filtered by
        drug or drug class.
    card_drug_map: passed straight through to Genotype.from_result_row,
        same as elsewhere in rules_engine.py.

    Returns a list of Genotype objects for this sample containing 
        nucleotide variant in multi-copy gene and gene copy number variants 
        added alongside the existing individual matches 
        (best match <= threshold)
    """
    if not geno_objs:
        return geno_objs

    organism = geno_objs[0].organism
    result = list(geno_objs)

    # intialise list of groups of Genotype objects, grouped by marker 
    # (for nucleotide variant in multi-copy gene) or gene (for gene copy number variant)
    nucl_groups = defaultdict(list)
    gene_groups = defaultdict(list)
    for g in geno_objs:
        if g.has_rule and not g.duplicated_row and g.rule.get('variation type') == 'Nucleotide variant detected in multi-copy gene':
            nucl_groups[g.marker_amrrules].append(g)
        if g.has_rule and not g.duplicated_row and g.rule.get('variation type') == 'Gene presence detected':
                gene_groups[g.gene_symbol].append(g)

    # first assess nucleotide variants
    if nucl_groups:
        base_nucl_rules = [
            r for r in rules
            if r.get('organism') == organism
            and r.get('variation type') == 'Nucleotide variant detected in multi-copy gene'
        ]
        parsed_nucl_rules = []
        for r in base_nucl_rules:
            marker, threshold = parse_multicopy_rule_mutation(r.get('mutation'), get_marker=True, gene=r.get('gene'))
            parsed_nucl_rules.append((marker, threshold, r))

        for marker, genos in nucl_groups.items():
            observed_copies = len(genos)
            possible_rules = [(threshold, r) for m, threshold, r in parsed_nucl_rules if m == marker]
            best_rule = _best_tier_rule(possible_rules, observed_copies)
            if best_rule:
                # duplicate the genotype object for the best rule, and add it to the result list
                new_geno = Genotype.from_result_row(genos[0], card_map=card_drug_map, rule=best_rule)
                # mark the row as a copy number row, so when we print to interpreted output
                # we can remove all the AMRFP info and just present the rule information
                new_geno.copy_number_row = True
                # store the original marker so we can match on it later
                new_geno.original_amrrules_marker = new_geno.marker_amrrules
                # update the marker to be a list of the markers, that is the length of the gene copies detected
                count = 1
                new_geno.marker_amrrules = marker
                while count < observed_copies:
                    count += 1
                    new_geno.marker_amrrules += f';{marker}'
                result.append(new_geno)

    # now look for gene copy number variants
    if gene_groups:
        base_gene_copy_rules = [
            r for r in rules
            if r.get('organism') == organism
            and r.get('variation type') == 'Gene copy number variant detected'
        ]
        parsed_gene_copy_rules = []
        for r in base_gene_copy_rules:
            marker, threshold = parse_multicopy_rule_mutation(r.get('mutation'), get_marker=True, gene=r.get('gene'))
            parsed_gene_copy_rules.append((marker, threshold, r))

        for gene, members in gene_groups.items():
            observed_copies = len(members)
            possible_rules = [(threshold, r) for m, threshold, r in parsed_gene_copy_rules if m == gene]
            best_rule = _best_tier_rule(possible_rules, observed_copies)
            if best_rule:
                new_geno = Genotype.from_result_row(members[0], card_map=card_drug_map, rule=best_rule)
                # members[0].variation_type was 'Gene presence detected' (copied
                # from the original item) - this new object represents a different
                # variation type entirely
                new_geno.variation_type = 'Gene copy number variant detected'
                # mark the row as a copy number row, so when we print to interpreted output
                # we can remove all the AMRFP info and just present the rule information
                new_geno.copy_number_row = True
                # store the original marker so we can match on it later
                new_geno.original_amrrules_marker = new_geno.marker_amrrules
                # update the marker to be a list of the markers, that is the length of the gene copies detected
                count = 1
                new_geno.marker_amrrules = gene
                while count < observed_copies:
                    count += 1
                    new_geno.marker_amrrules += f';{gene}'
                result.append(new_geno)

    return result

def apply_combination_rules(geno_objs, rules, card_drug_map):
    if not geno_objs:
        return geno_objs
    organism = geno_objs[0].organism
    result = list(geno_objs)

    for g in geno_objs:
        # get all the solo rule ids to compare against
        solo_rule_ids = [g.ruleID for g in geno_objs if g.has_rule and not g.duplicated_row]
        if not solo_rule_ids:
            continue
    # extract all the combination rules for this organism
    combo_rules = get_combination_rules(rules, organism)
    for rule in combo_rules:
        ruleID_logic = rule.get('gene')
        # for each combo rule, check if we have all the solo rules required to satisfy the logic string
        if evaluate_logic_string(ruleID_logic, solo_rule_ids):
            rules_in_logic = set(re.findall(r'\b\w+\b', ruleID_logic))
            # extract the relevant rules
            matching_objs = [g for g in geno_objs if g.ruleID in rules_in_logic]

            # build the combo's marker string by substituting each ruleID in the logic string with 
            # its marker, preserving the logic structure
            combo_marker = ruleID_logic
            for g in matching_objs:
                combo_marker = re.sub(rf'\b{re.escape(g.ruleID)}\b', g.marker_amrrules, combo_marker)
            # now we need to create a new genotype object for this combination rule
            new_geno = Genotype.from_result_row(matching_objs[0], card_map=card_drug_map, rule=rule)
            new_geno.combo_rule_row = True
            new_geno.marker_amrrules = combo_marker
            # store the list of ruleIDs that make up the combination rule
            new_geno.combo_rule_components = rules_in_logic
            # add it to the list of results to return
            result.append(new_geno)
            # update the subcomponent parts to indicate they are part of a combo rule
            #for g in matching_objs:
            #    result.remove(g)

    return result