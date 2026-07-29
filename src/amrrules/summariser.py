from amrrules.resources import ResourceManager as rm
from amrrules.utils import CATEGORY_ORDER, PHENOTYPE_ORDER, EVIDENCE_GRADE_ORDER
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
    
    def summarise_rules(self, no_rule_interpretation, combo_rules, class_summary=None):
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
            # move the partial call to the ruleID col and set drug and class to '-' to avoid confusion
            if self.drug_class == 'partial':
                self.ruleIDs = 'none (partial hits)'
                self.drug_class = '-'
                self.drug = '-'
            return
        if self.drug_class == 'antibiotic efflux':
            self.drug = '(n/a)'
        
        # otherwise, continue on

        # first, grab all the individual ruleIDs that have been applied to this drug or drug class
        solo_rule_ids = [g.ruleID for g in geno_objs
                    if getattr(g, "ruleID", None) not in (None, "-")]
        # Set the rule IDs in the output, or '-' if none were found
        self.ruleIDs = ";".join(sorted(solo_rule_ids)) if solo_rule_ids else "-"

        # If we have combination rules, we need to evaluate them to see if any apply.
        # But this should only be evaluated if we have rules that are being applied
        rules_overriden_by_combo = set()
        combo_rule_id_matches = []
        rules_to_assess = []
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
                    rules_overriden_by_combo.update(rules_in_logic)
                    # add to the list of applied combo rules for printing to output
                    combo_rule_id_matches.append(rule.get('ruleID'))
                    # add this rule to the list of rules to assess for interpretation
                    rules_to_assess.append(rule)

        # Set combo rule IDs in the output, or '-' if none were found
        if len(combo_rule_id_matches) == 0:
            self.combo_rules = '-'
        else:
            self.combo_rules = ";".join(combo_rule_id_matches)

        # Update our rules to assess by only including solo individual rules that were not overridden by a combo rule
        for g in geno_objs:
            if g.ruleID not in rules_overriden_by_combo and g.ruleID not in (None, "-"):
                rules_to_assess.append(g.rule)

        # First, set overall WT/NWT status based on the rules we need to assess
        phenotypes = [r['phenotype'] for r in rules_to_assess if 'phenotype' in r]
        self.phenotype = self._get_max_value(phenotypes, PHENOTYPE_ORDER)

        # Determine the overall clinical category based on the rules we need to assess
        clinical_categories = [r['clinical category'] for r in rules_to_assess if 'clinical category' in r]
        self.category = self._get_max_value(clinical_categories, CATEGORY_ORDER)

        # Finally, set the overall evidence grade. This is highest evidence grade linked to any rules matching our highest clinical category.
        evidence_grades = [r['evidence grade'] for r in rules_to_assess if 'evidence grade' in r and r['clinical category'] == self.category]
        self.evidence_grade = self._get_max_value(evidence_grades, EVIDENCE_GRADE_ORDER)

        # alright, but depending on our no_rule_interpretation setting, we may need to override the category and phenotype values
        # only matters if we have markers with no rules
        if self.markers_with_norule != '-':
            if no_rule_interpretation == 'none' or no_rule_interpretation == 'nwt':
                # for the category, if we have any nwt markers without rules, then we can't interpret
                # what this means in combination with an S marker, so set to '-'
                # however if the rule says 'R', then we can keep the R
                if self.category == 'S':
                    self.category = '-'
                # evidence grade also gets switched to 'none' regardless
                self.evidence_grade = 'none'
                # we change the phenotype based on whether its none or nwt
                if no_rule_interpretation == 'none':
                # for the phenotype, if we have any nwt markers without rules, then we can't interpret, so set to '-'
                    self.phenotype = '-'
                elif no_rule_interpretation == 'nwt':
                    # in this case, if our rule markers state we have a wt phenotype, but we have nwt markers with no rule
                    # we override the penotype to be nwt
                    # if we had markers with rules that were nwt R, we stay nwt anyway
                    self.phenotype = 'nonwildtype'
            if no_rule_interpretation == 'nwtS':
                # if we have any nwt markers without rules, we're calling nwt and S
                # but we want the evidence grade to be 'none' to reflect the fact
                # that the call is being made using markers with no rules
                self.evidence_grade = 'none'

        # if efflux, then set clinical category to '-'
        if self.drug_class == 'antibiotic efflux':
            self.category = '-'

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

    # Helper to get max by order list
    @staticmethod
    def _get_max_value(values, order):
        valid_values = [v for v in values if v in order]
        if not valid_values:
            return None
        return max(valid_values, key=lambda v: order.index(v))


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
                # extract combo rules for this drug_class
                # to get the list of possible combo rules to evaluate, we need to extract all 'Combination' rules for this organism
                combo_rules = [r for r in rules if r.get('organism') == summary_entry.organism and r.get('variation type') == 'Combination']
                # then need to further filter to include only combo rules that apply to the drug class we're assessing
                combo_rules = [r for r in combo_rules if summary_entry.drug_class in r.get('drug class', '')]
                # determine the highest category/pheno/evidence grade for this drug_class
                summary_entry.summarise_rules(no_rule_interpretation, combo_rules)
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
                    # but take into account the rules for the drug class
                    combo_rules = [r for r in rules if r.get('organism') == summary_entry.organism and r.get('variation type') == 'Combination']
                    # then need to further filter to include only combo rules that apply to either the drug or class we're assessing
                    combo_rules = [r for r in combo_rules if summary_entry.drug in r.get('drug', '') or summary_entry.drug_class in r.get('drug class', '')]
                    summary_entry.summarise_rules(no_rule_interpretation, combo_rules, class_summary=master_class_entry)
                    # add it to our list
                    summary_entry_list.append(summary_entry)
            summary_entry_dict[sample_name] = order_summary_objs(summary_entry_list)
    
    return(summary_entry_dict)