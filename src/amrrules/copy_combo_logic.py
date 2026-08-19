from collections import defaultdict
import ast
import re
from amrrules.rules_io import parse_multicopy_rule_mutation, get_combination_rules, evaluate_logic_string
from amrrules.genotype_parser import Genotype


def _simplify_logic_expression(node, marker_by_ruleid):
    """
    Recursively simplify a combo rule's ruleID_logic string, parsed with 
    `and`/`or` in place of &/|) down to a marker string, keeping only 
    branches whose ruleID was actually detected. 
    Returns None if this node has no detected component at all.
    """
    if isinstance(node, ast.Name):
        return marker_by_ruleid.get(node.id)

    if isinstance(node, ast.BoolOp):
        simplified_children = [_simplify_logic_expression(v, marker_by_ruleid) for v in node.values]
        present = [c for c in simplified_children if c is not None]

        if isinstance(node.op, ast.And):
            # AND requires every operand to be present for the overall
            # expression to have evaluated True - if one's missing, this
            # rule shouldn't have matched at all
            if len(present) != len(simplified_children):
                raise ValueError("AND branch missing a detected marker - rule shouldn't have matched")
            return " & ".join(present)

        if isinstance(node.op, ast.Or):
            if not present:
                return None
            if len(present) == 1:
                return present[0]
            return "(" + " | ".join(present) + ")"

    raise ValueError(f"Unsupported logic expression node: {ast.dump(node)}")


def simplify_combo_marker(ruleID_logic, matching_objs):
    """
    Given a combo rule's ruleID_logic string (e.g. 'A & B & (C | D) & E')
    and the Genotype objects actually detected/matched, return a marker
    string with undetected OR-branches dropped entirely - e.g.
    'A_marker & B_marker & C_marker & E_marker' if D wasn't detected.
    """
    marker_by_ruleid = {g.ruleID: g.marker_amrrules for g in matching_objs}
    python_logic = ruleID_logic.replace('&', ' and ').replace('|', ' or ')
    tree = ast.parse(python_logic, mode='eval').body
    result = _simplify_logic_expression(tree, marker_by_ruleid)
    if result is None:
        raise ValueError(f"No detected markers satisfy logic string: {ruleID_logic}")
    return result


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

    # get all the solo rule ids to compare against
    solo_rule_ids = [g.ruleID for g in geno_objs if g.has_rule and not g.duplicated_row]
    if not solo_rule_ids:
        return result
    # extract all the combination rules for this organism
    combo_rules = get_combination_rules(rules, organism)
    for rule in combo_rules:
        ruleID_logic = rule.get('gene')
        # for each combo rule, check if we have all the solo rules required to satisfy the logic string
        if evaluate_logic_string(ruleID_logic, solo_rule_ids):
            rules_in_logic = set(re.findall(r'\b\w+\b', ruleID_logic))
            # extract the relevant rules
            matching_objs = [g for g in geno_objs if g.ruleID in rules_in_logic]
            # build the combo's marker string, dropping any undetected
            # OR-branches entirely rather than leaving their ruleIDs in place
            combo_marker = simplify_combo_marker(ruleID_logic, matching_objs)
            # now we need to create a new genotype object for this combination rule
            new_geno = Genotype.from_result_row(matching_objs[0], card_map=card_drug_map, rule=rule)
            new_geno.combo_rule_row = True
            new_geno.marker_amrrules = combo_marker
            # store the list of ruleIDs that make up the combination rule
            new_geno.combo_rule_components = rules_in_logic
            # add it to the list of results to return
            result.append(new_geno)

    return result