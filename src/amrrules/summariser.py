from amrrules.resources import ResourceManager as rm
from amrrules.utils import PAIRWISE_TABLE, DEFAULT_COMBINE_TABLE, IMPOSSIBLE, ImpossibleCombination, _normalize_call, EVIDENCE_GRADE_ORDER
from amrrules.rules_io import parse_multicopy_rule_mutation
from collections import defaultdict
import re

class SummaryEntry:
    """
    Class to summarise genotype rows by drug or drug class for a given sample.
    """

    def __init__(self, sample_name, genotype_objects):
        
        # these are all the columns that are going to be in the output
        # these values will apply for a particular drug or drug class
        self.sample_name = sample_name
        self.geno_objs = genotype_objects
        self.drug = genotype_objects[0].drug
        self.drug_class = genotype_objects[0].drug_class
        self.organism = genotype_objects[0].organism # this will be the same for all objects

        # if we're just working with a drug_class, then drug should be set to (all) to make clear
        # that this applies to the whole class
        if self.drug == '-' and self.drug_class != '-' and self.drug_class != 'unassigned markers':
            self.drug = '(all)'

        # these are also columns that we're going to set with the below functions
        self.category = None
        self.gene_context = None
        self.phenotype = None
        self.evidence_grade = None
        self.markers_rule_nonS = None
        self.markers_with_norule = None
        self.markers_s = None
        self.ruleIDs = None
        self.combo_rules = None
    
    def summarise_rules(self, no_rule_interpretation, combo_rules, multi_copy_rules, class_summary=None):
        """Compute summary values based on geno_objs."""

        # full object list creation, if we've got a drug class of objs also to consider
        if class_summary:
            geno_objs = self.geno_objs + class_summary.geno_objs
        else:
            geno_objs = self.geno_objs

        # if our class is 'unassigned markers' or 'partial', then we have no category/phenotype/evidence
        # so just set these values and exit
        if self.drug_class in ['unassigned markers', 'partial']:
            self.category = '-'
            self.phenotype = '-'
            self.evidence_grade = '-'
            self.ruleIDs = '-'
            self.combo_rules = '-'
            # move the partial call to the ruleID col and set drug and class to '-' to avoid confusion
            if self.drug_class == 'partial':
                self.ruleIDs = 'none (partial hits)'
                self.drug_class = '-'
                self.drug = '-'
            return

        # if we have no rules to apply, then we can't interpret
        # so set the values to match the no_rule_interpretation setting, and exit
        if not any(g.has_rule for g in geno_objs):
            # no rules to apply, therefore these values are '-'
            self.ruleIDs = '-'
            self.combo_rules = '-'
            if no_rule_interpretation == 'none':
                self.category = '-'
                self.phenotype = '-'
                self.evidence_grade = 'none'
            elif no_rule_interpretation == 'nwt':
                self.category = '-'
                self.phenotype = 'nonwildtype'
                self.evidence_grade = 'none'
            elif no_rule_interpretation == 'nwtS':
                self.category = 'S'
                self.phenotype = 'nonwildtype'
                self.evidence_grade = 'none'
            elif no_rule_interpretation == 'nwtR':
                self.category = 'R'
                self.phenotype = 'nonwildtype'
                self.evidence_grade = 'none'
            if self.drug_class == 'antibiotic efflux':
                #override as we can't say anything meaningful for efflux
                self.category = '-'
                self.phenotype = '-'
                self.evidence_grade = '-'
                self.drug = '(n/a)'
            return
        
        # otherwise, continue on
        # first, grab all the individual ruleIDs that have been applied to this drug or drug class
        solo_rule_ids = set(g.ruleID for g in geno_objs
                    if getattr(g, "ruleID", None) not in (None, "-"))
        # use this to keep track of any rules that are going to be overridden by a combination or multi-copy gene rule
        rules_to_be_overriden = set()
        # use this to keep track of rules that need to be assessed for the final interpretation. This will be a combination of solo rules, combination rules, and any multi-copy rules that apply
        rules_to_assess = []

        # now, check to see if we have any multi-copy rules for this organism + drug/drug class combination
        if multi_copy_rules:
            # if we have multi copy rules for this drug, check to see if any of our geno_objects would match any of the multi-copy rules
            # store the geno objects together in a dict by their marker_amrrules value, so we can check what the total length of this dict is
            nucl_variant_multicopy_genos = {}
            for g in geno_objs:
                # first lets look at geno objects that have rules, and are multi-copy nucleotide variants
                if g.has_rule and g.rule.get('variation type') == 'Nucleotide variant detected in multi-copy gene':
                    nucl_variant_multicopy_genos.setdefault(g.marker_amrrules, []).append(g)
            #now go through each marker in the dict, and count the number of objects that exist for that marker
            for marker, marker_genos in nucl_variant_multicopy_genos.items():
                # this is the total number of copies we've observed
                observed_copy_count = len(marker_genos)
                # now go through the multi-copy rules and find any that match this marker and have a threshold <= observed_copy_count
                matching_rules = [r for r in multi_copy_rules if r.get('marker_amrrules') == marker and r.get('threshold') <= observed_copy_count]
                # if we have any matching rules, then we need to find the one with the highest threshold
                if matching_rules:
                    # sort the matching rules by threshold, and take the last one (highest threshold)
                    best_rule = sorted(matching_rules, key=lambda r: r.get('threshold'))[-1]
                    # this is now the best rule for this set of objects. We want to add this rule to the list of ruleIDs for us to assess
                    rules_to_assess.append(best_rule)
                    # add the individual rules for these markers to the list of rules to be overridden, as the multi-copy rule overrides the individual rules
                    for g in marker_genos:
                        rules_to_be_overriden.add(g.ruleID)

        # Set the rule IDs in the output, or '-' if none were found
        self.ruleIDs = ";".join(sorted(solo_rule_ids)) if solo_rule_ids else "-"

        # If we have combination rules, we need to evaluate them to see if any apply.
        combo_rule_id_matches = set()
        if solo_rule_ids and combo_rules:
            for rule in combo_rules:
                ruleID_logic = rule.get('gene')
                # ruleID logic is a string of the form "gene1 & gene2 | gene3"
                # we need to replace the & with a python 'and' and the | with a python 'or'
                matched_combo = self._evaluate_logic_string(ruleID_logic, solo_rule_ids)
                if matched_combo:
                    # extract the individual rule IDs so we can exclude these rules from our
                    #interpretation logic later, as the combo rule overrides the individual rules
                    rules_in_logic = set(re.findall(r'\b\w+\b', ruleID_logic))
                    rules_to_be_overriden.update(rules_in_logic)
                    # add to the list of applied combo rules for printing to output
                    combo_rule_id_matches.add(rule.get('ruleID'))
                    # add this rule to the list of rules to assess for interpretation
                    rules_to_assess.append(rule)

        # Set combo rule IDs in the output, or '-' if none were found
        if len(combo_rule_id_matches) == 0:
            self.combo_rules = '-'
        else:
            self.combo_rules = ";".join(combo_rule_id_matches)

        # Update our rules to assess by only including solo individual rules that were not overridden by a combo rule
        for g in geno_objs:
            if g.ruleID not in rules_to_be_overriden and g.ruleID not in (None, "-"):
                rules_to_assess.append(g.rule)

        # extract all the calls for the rules we need to assess
        calls = [(r['phenotype'], r['clinical category']) for r in rules_to_assess]
        # determine the overall call for this set of rules
        rule_call = self.combine_many(calls)

        # if we have markers with no rule, then update the rule call based on our default setting
        if self.markers_with_norule != '-':
            rule_call = self.apply_norule_default(rule_call, no_rule_interpretation)

        self.phenotype, self.category = rule_call

        # now get the evidence grade for the final call, based only on the rules
        # that match our final clinical category
        evidence_grades = [
        r['evidence grade'] for r in rules_to_assess
        if 'evidence grade' in r and r.get('clinical category') == self.category]

        # if there are no evidence grades, then we had no rules that match the final category
        if not evidence_grades:
            self.evidence_grade = 'none'
        else:
            self.evidence_grade = max(evidence_grades, key=lambda v: EVIDENCE_GRADE_ORDER.index(v))

        if self.drug_class == 'antibiotic efflux':
                    #override as we can't say anything meaningful for efflux
                    self.category = '-'
                    self.phenotype = '-'
                    self.evidence_grade = '-'
                    self.drug = '(n/a)'

        return

    def combine_calls(self, call1, call2):
        """Look up two (phenotype, category) calls in PAIRWISE_TABLE."""
        call1 = _normalize_call(call1)
        call2 = _normalize_call(call2)
        result = PAIRWISE_TABLE[call1][call2]
        if result is IMPOSSIBLE:
            raise ImpossibleCombination(f"{call1} + {call2} is not a valid combination")
        return result
 
    def combine_many(self, calls):
        """
        Combine a list of (phenotype, category) calls into one, folding all
        category == 'S' calls together first, then folding the remaining
        (non-S) calls in one at a time against the running result.
    
        This ordering should never hit an *np cell by design. If it does, a
        warning is printed and that call is skipped (running result kept
        unchanged) rather than crashing - so the source rules can be fixed.
        """
        if not calls:
            return None
        calls = [_normalize_call(c) for c in calls]
    
        s_calls = [c for c in calls if c[1] == 'S']
        non_s_calls = [c for c in calls if c[1] != 'S']
        ordered_calls = s_calls + non_s_calls
    
        result = ordered_calls[0]
        for call in ordered_calls[1:]:
            try:
                result = self.combine_calls(result, call)
            except ImpossibleCombination as e:
                print(
                    f"WARNING: impossible combination hit while combining rules "
                    f"({e}). This should not happen by design - check the source "
                    f"rules. Skipping this call and keeping the current result."
                )
        return result

    def apply_norule_default(self, rule_call, no_rule_interpretation):
        """
        Combine the winning rule-based call with the no_rule_interpretation
        default, but ONLY if the sample actually has markers with no matching
        rule. If there are none, the rule-based call passes through unchanged.
        """
        
        default_row = DEFAULT_COMBINE_TABLE[no_rule_interpretation]
        return default_row[_normalize_call(rule_call)]

    def set_markers(self, flag_core, class_summary=None):
        
        # for each object, extract the marker and place it into the correct
        # list based on whether it has a rule, no rule, or is wildtype
        markers_rule_nonS = []
        markers_with_norule = []
        markers_s = []

        # first loop through the markers for the drug
        for g in self.geno_objs:
            # set the marker to be the amrrules formatted version
            marker = g.marker_amrrules
            if g.has_rule:
                # only label core genes if it's a core context with gene presence variation type
                if g.gene_context == 'core' and g.variation_type == 'Gene presence detected' and flag_core:
                    marker = marker + " (core)"
                if g.clinical_category == 'S':
                    markers_s.append(marker)
                else:
                    markers_rule_nonS.append(marker)
            else:
                markers_with_norule.append(marker)
        
        # now loop through the markers for the class, if there is one
        if class_summary:
            for g in class_summary.geno_objs:
                # set the marker to be the amrrules formatted version
                marker = g.marker_amrrules
                if g.has_rule:
                    # only label core genes if it's a core context with gene presence variation type
                    if g.gene_context == 'core' and g.variation_type == 'Gene presence detected' and flag_core:
                        marker = marker + " (core)"
                    # append marker only if it's not already present
                    if g.clinical_category == 'S' and marker not in markers_s:
                        markers_s.append(marker)
                    elif g.clinical_category != 'S' and marker not in markers_rule_nonS:
                        markers_rule_nonS.append(marker)
                # append marker only if it's not already present
                elif marker not in markers_with_norule:
                    markers_with_norule.append(marker)

        self.markers_rule_nonS = ';'.join(markers_rule_nonS) or '-'
        self.markers_with_norule = ';'.join(markers_with_norule) or '-'
        self.markers_S = ';'.join(markers_s) or '-'

    @staticmethod
    def _evaluate_logic_string(logic_string, id_list):
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

def order_summary_objs(objs):
    """
    Sort a list of summaryEntry objects first by drug_class (alphabetically, with 'antibiotic efflux', followed by 'unassigned markers', with 'partial' last),
    then by drug (alphabetically, '-' last).
    """
    def drug_class_sort_key(obj):
        drug_class_lower = getattr(obj, "drug_class", "").lower()
        rule_ids = getattr(obj, "ruleIDs", "")
        if rule_ids == "none (partial hits)":
            return (3, "")  # Last
        elif drug_class_lower == "unassigned markers":
            return (2, "")  # Second last
        elif drug_class_lower == "antibiotic efflux":
            return (1, "")  # Third last
        else:
            return (0, drug_class_lower)  # Alphabetical for all others
    
    sorted_list = sorted(
        objs,
        key=lambda o: (
            drug_class_sort_key(o),
            # For drug: alphabetical, with '-' last
            (getattr(o, "drug", "").lower() == "-", getattr(o, "drug", "").lower())
        )
    )

    return sorted_list

def get_combination_rules(rules, organism, drug_class, drug=None):
    """
    Extracts combination rules or multi-copy rules for a given organism and drug/drug class.
    Always filters first for drug class as this will always be provided.
    Drug then needs to be added on if a specific drug is provided, in case there are specific rules for the drug

    Returns a tuple of (combo_rules, multi_copy_rules) where each is a list of rules that match the organism and drug/drug class.
    """
    base_rules = [
        r for r in rules
        if r.get('organism') == organism
        and r.get('variation type') in (
            'Combination',
            'Nucleotide variant detected in multi-copy gene',
            'Gene copy number variant detected',
        )
    ]
    matching_rules = [r for r in base_rules if r.get('drug class') == drug_class]
    if drug is not None:
        matching_rules.extend([r for r in base_rules if drug in r.get('drug', '')])

    combo_rules = [r for r in matching_rules if r.get('variation type') == 'Combination']
    multi_copy_rules = [r for r in matching_rules if r.get('variation type') != 'Combination']

    # for the multi-copy rules, we need to parse each rule's mutation to extract the base mutation and threshold, and add these as new keys in the rule dict
    for rule in multi_copy_rules:
        mutation = rule.get('mutation')
        if mutation:
            marker, threshold = parse_multicopy_rule_mutation(mutation, get_marker=True, gene=rule.get('gene'))
            rule['threshold'] = threshold
            # this is the key we're going to use to match
            rule['marker_amrrules'] = marker

    return combo_rules, multi_copy_rules

def create_summary_dict(grouped_by_sample, rules, flag_core, no_rule_interpretation):

    summary_entry_dict = {} # key: sample name, value: list of summary entry objs
    for sample_name, genotypes in grouped_by_sample.items():
        summary_entry_list = []
        # for the sample, we need to group by drug class, and then by drug within that

        sample_groups = defaultdict(lambda: defaultdict(list))
        for g in genotypes:
            sample_groups[g.drug_class][g.drug].append(g)
        for drug_class in sample_groups.keys():
            # for each drug_class, we first need to apply a summary entry at the class level
            # if the class level exists
            class_level_hits = sample_groups[drug_class].get('-', None)
            # initialise our master class entry as None, will be filled later
            master_class_entry = None
            if class_level_hits:
                summary_entry = SummaryEntry(sample_name, class_level_hits)
                # assign markers with, without rules, and wt markers
                summary_entry.set_markers(flag_core)
                combo_rules, multi_copy_rules = get_combination_rules(rules, summary_entry.organism, summary_entry.drug_class)
                summary_entry.summarise_rules(no_rule_interpretation, combo_rules, multi_copy_rules)
                # this is our master entry for this drug_class, so save it
                master_class_entry = summary_entry
                # add it to our list
                summary_entry_list.append(summary_entry)
            
            # otherwise now we're in a specific drug for the class
            # we need to make sure that the interpretation of this drug doesn't conflict with the class level rules
            for drug in sample_groups[drug_class].keys():
                if drug != '-':
                    # create our summary entry
                    summary_entry = SummaryEntry(sample_name, sample_groups[drug_class][drug])
                    # before assigning markers to columns
                    # we want to remove any duplicated row markers from the class level
                    # assign markers
                    summary_entry.set_markers(flag_core, class_summary=master_class_entry)
                    # determine highest category/pheno/evidence grade for this drug, including combo rules (if any)
                    # but take into account any combination or multi copy rules for the drug class or drug
                    combo_rules, multi_copy_rules = get_combination_rules(rules, summary_entry.organism, summary_entry.drug_class, summary_entry.drug)
                    summary_entry.summarise_rules(no_rule_interpretation, combo_rules, multi_copy_rules, class_summary=master_class_entry)
                    # add it to our list
                    summary_entry_list.append(summary_entry)
            summary_entry_dict[sample_name] = order_summary_objs(summary_entry_list)
    
    return(summary_entry_dict)