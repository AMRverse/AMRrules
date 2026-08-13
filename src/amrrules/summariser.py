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
    
    #def summarise_rules(self, no_rule_interpretation, combo_rules, class_summary=None, flag_core=False):
    def summarise_rules(self, no_rule_interpretation, class_summary=None, flag_core=False):

        """Compute summary values based on geno_objs."""

        # full object list creation, if we've got a drug class of objs also to consider
        if class_summary:
            geno_objs = []
            # go through the class summary geno objects
            # only add the ones that have different markers to the drug level objects
            for g in class_summary.geno_objs:
                if g.marker_amrrules not in {x.marker_amrrules for x in self.geno_objs}:
                    geno_objs.append(g)
            # now add the drug level objects
            geno_objs.extend(self.geno_objs)
        else:
            # otherwise our objects to parse are simply the drug level objects
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
            self.set_markers(geno_objs, {}, flag_core=flag_core)
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
            self.set_markers(geno_objs, {}, flag_core=flag_core)
            return
        
        # otherwise, continue on
        # first, grab all the individual ruleIDs that have been applied to this drug or drug class
        solo_rule_ids = set()
        combo_rule_ids = set()
        for g in geno_objs:
            if getattr(g, "ruleID", None) in (None, "-"):
                continue
            if getattr(g, "combo_rule_row", False):
                combo_rule_ids.add(g.ruleID)
            else:
                solo_rule_ids.add(g.ruleID)

        self.combo_rules = ";".join(sorted(combo_rule_ids)) if combo_rule_ids else '-'
        # use this to keep track of any rules that are going to be overridden by a combination or multi-copy gene rule
        rules_to_be_overriden = set()
        # use this to keep track of rules that need to be assessed for the final interpretation. This will be a combination of solo rules, combination rules, and any multi-copy rules that apply
        rules_to_assess = []

        # a multi-copy gene object overrides any individual objects that are part of the multi-copy rule
        marker_for_objs_to_remove = None
        copy_number_override = False
        combo_markers_to_add = defaultdict(list)
        for g in geno_objs:
            if g.copy_number_row:
                copy_number_override = True
                # this is the master row, so remove any other geno objects that have the same marker
                marker_for_objs_to_remove = g.original_amrrules_marker
            if getattr(g, "combo_rule_row", False):
                # this is a combination rule, so add the individual rules that make up the combination
                # to the list of rules to be overriden (excluded) when generating the final call
                # only do this if our rules are already present in the solo_rule_ids list, as there may be instances where a combination rule is applied, but the individual rules that make up the combo are not present in our drug specific list, as the individual components are relevant to a different drug
                if set(g.combo_rule_components).issubset(solo_rule_ids):
                    rules_to_be_overriden.update(g.combo_rule_components)
                # if this is the case, then not only do we want to add the rule IDs, we want 
                # to add the individual markers to a dict of markers to add to the final marker strings
                else:
                    solo_rule_ids.update(g.combo_rule_components)
                    combo_markers_to_add.setdefault(g.clinical_category, set()).update(re.findall(r"[^&|()\s]+", g.marker_amrrules))

        # now remove any geno objects that have the same marker as the multi-copy row, 
        # and are not the multi-copy row itself
        to_keep = []
        for g in geno_objs:
            if copy_number_override and (g.marker_amrrules == marker_for_objs_to_remove and not g.copy_number_row):
                rules_to_be_overriden.add(g.ruleID)
                solo_rule_ids.discard(g.ruleID)
            else:
                to_keep.append(g)
        geno_objs = to_keep


        # If we have combination rules, we need to evaluate them to see if any apply.
        """
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
        """

        # now set all the markers
        # only assess the genotype objects that aren't being overridden
        self.set_markers(geno_objs, combo_markers_to_add, flag_core=flag_core)

        # Set the rule IDs in the output, or '-' if none were found
        # doing this here so we exclude any rule IDs from markers that have been collapsed
        # into a multi copy rule (but keeping individual rule IDs that make up combination rules)
        self.ruleIDs = ";".join(sorted(solo_rule_ids)) if solo_rule_ids else "-"

        # update the rules to assess list to remove any rules that are 
        # being overridden by combination or multi-copy rules
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

    def set_markers(self, geno_objs, combo_markers_to_add, flag_core):
        
        # for each object, extract the marker and place it into the correct
        # list based on whether it has a rule, no rule, or is wildtype
        markers_rule_nonS = []
        markers_with_norule = []
        markers_s = []

        # set up a list of processed markers,
        # so we can exclude any duplicate markers from
        # the class level, that are being inherited
        processed_markers = []

        # first loop through the markers for the drug
        for g in geno_objs:
            # set the marker to be the amrrules formatted version
            marker = g.marker_amrrules
            processed_markers.append(marker)
            if g.has_rule:
                # skip markers that are part of combo rules, as the individual
                # markers will be included in the correct marker string already
                if getattr(g, "combo_rule_row", False):
                    continue
                # only label core genes if it's a core context with gene presence variation type
                if g.gene_context == 'core' and g.variation_type == 'Gene presence detected' and flag_core:
                    marker = marker + " (core)"
                if g.clinical_category == 'S':
                    markers_s.append(marker)
                else:
                    markers_rule_nonS.append(marker)
            else:
                markers_with_norule.append(marker)
        
        # now loop through any combo markers that need to be added
        for category, markers in combo_markers_to_add.items():
            for marker in markers:
                if marker in processed_markers:
                    continue
                processed_markers.append(marker)
                if category == 'S':
                    markers_s.append(marker)
                else:
                    markers_rule_nonS.append(marker)

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
    Extracts combination rules for a given organism and drug/drug class.
    Always filters first for drug class as this will always be provided.
    Drug then needs to be added on if a specific drug is provided, in case there are specific rules for the drug

    Returns a list of combination rules that match the organism and drug/drug class.
    """
    base_rules = [
        r for r in rules
        if r.get('organism') == organism
        and r.get('variation type') in ('Combination')
    ]
    matching_rules = [r for r in base_rules if r.get('drug class') == drug_class]
    if drug is not None:
        matching_rules.extend([r for r in base_rules if drug in r.get('drug', '')])

    combo_rules = [r for r in matching_rules if r.get('variation type') == 'Combination']


    return combo_rules

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
                #combo_rules = get_combination_rules(rules, summary_entry.organism, summary_entry.drug_class)
                #summary_entry.summarise_rules(no_rule_interpretation, combo_rules, flag_core=flag_core)
                summary_entry.summarise_rules(no_rule_interpretation, flag_core=flag_core)
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
                    # determine highest category/pheno/evidence grade for this drug, including combo rules (if any)
                    # but take into account any combination rules for the drug class or drug
                    #combo_rules = get_combination_rules(rules, summary_entry.organism, summary_entry.drug_class, summary_entry.drug)
                    #summary_entry.summarise_rules(no_rule_interpretation, combo_rules, class_summary=master_class_entry, flag_core=flag_core)
                    summary_entry.summarise_rules(no_rule_interpretation, class_summary=master_class_entry, flag_core=flag_core)
                    # add it to our list
                    summary_entry_list.append(summary_entry)
            summary_entry_dict[sample_name] = order_summary_objs(summary_entry_list)
    
    return(summary_entry_dict)