from amrrules.resources import ResourceManager as rm
from amrrules.utils import CATEGORY_ORDER, PHENOTYPE_ORDER, EVIDENCE_GRADE_ORDER
from collections import defaultdict

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
    
    def summarise_rules(self, class_summary=None):
        """Compute summary values based on geno_objs."""

        # Helper to get max by order list
        def get_max_value(values, order):
            valid_values = [v for v in values if v in order]
            if not valid_values:
                return None
            return max(valid_values, key=lambda v: order.index(v))

        # full object list creation, if we've got a drug class of objs also to consider
        if class_summary:
            self.geno_objs = self.geno_objs + class_summary.geno_objs
        # Extract values from genotype objects
        #categories = [g.clinical_category for g in self.geno_objs if hasattr(g, 'clinical_category')]
        phenotypes = [g.phenotype for g in self.geno_objs if hasattr(g, 'phenotype')]
        #evidence_grades = [g.evidence_grade for g in self.geno_objs if hasattr(g, 'evidence_grade')]

        # update evidence grade to be linked to the evidence for the highest category call 
        # (eg if oqx is S with rule, but there is gyrA marker with no rule that's R, category is R but evidence grade is very low)
        best_obj = max(self.geno_objs, key=lambda o: (
            CATEGORY_ORDER.index(o.clinical_category),
            EVIDENCE_GRADE_ORDER.index(o.evidence_grade)
            ))

        # Set the “maximum” according to ordering
        self.category = best_obj.clinical_category
        self.phenotype = get_max_value(phenotypes, PHENOTYPE_ORDER)
        self.evidence_grade = best_obj.evidence_grade
    
    def set_ruleIDs_and_combo(self, combo_rules):
        # deal with the ruleIDs
        rule_ids = {g.ruleID for g in self.geno_objs
            if getattr(g, "ruleID", None) not in (None, "-")}
        # sets ruleIDs to '-' if there are none at all
        self.ruleIDs = ";".join(sorted(rule_ids)) if rule_ids else "-"

        matched_combo_rules = []
        # for each rule, we need to extract the ruleID logic and check if it matches our ruleIDs
        for rule in combo_rules:
            ruleID_logic = rule.get('gene')
            # ruleID logic is a string of the form "gene1 & gene2 | gene3"
            # we need to replace the & with a python 'and' and the | with a python 'or'
            matched_combo = self._evaluate_logic_string(ruleID_logic, rule_ids)
            if matched_combo:
                # add to the list of matched combos
                matched_combo_rules.append(ruleID_logic)
        # if we didn't find any matching combo rules, then we return None
        if len(matched_combo_rules) == 0:
            self.combo_rules = '-'
        else:
            self.combo_rules = ";".join(matched_combo_rules)
    
    def set_markers(self, flag_core):
        # for each object, extract the marker and place it into the correct
        # list based on whether it has a rule, no rule, or is wildtype
        markers_rule_nonS = []
        markers_with_norule = []
        markers_s = []
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
                #if g.clinical_category == 'S':
                #    markers_s.append(marker)
                #else:
                markers_with_norule.append(marker)

        self.markers_rule_nonS = ';'.join(markers_rule_nonS) or '-'
        self.markers_with_norule = ';'.join(markers_with_norule) or '-'
        self.markers_S = ';'.join(markers_s) or '-'

    def _evaluate_logic_string(logic_string, id_list):
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


def order_summary_objs(objs):
    """
    Sort a list of summaryEntry objects first by drug_class (alphabetically, with 'antibiotic efflux' second last and 'unassigned markers' last),
    then by drug (alphabetically, '-' last).
    """
    def drug_class_sort_key(drug_class):
        drug_class_lower = drug_class.lower()
        if drug_class_lower == "unassigned markers":
            return (2, "")  # Last
        elif drug_class_lower == "antibiotic efflux":
            return (1, "")  # Second last
        else:
            return (0, drug_class_lower)  # Alphabetical for all others
    
    sorted_list = sorted(
        objs,
        key=lambda o: (
            drug_class_sort_key(getattr(o, "drug_class", "")),
            # For drug: alphabetical, with '-' last
            (getattr(o, "drug", "").lower() == "-", getattr(o, "drug", "").lower())
        )
    )

    return sorted_list

def create_summary_dict(grouped_by_sample, rules, flag_core):

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
            # initialise our claster class entry as None, will be filled later
            master_class_entry = None
            if class_level_hits:
                summary_entry = SummaryEntry(sample_name, class_level_hits)
                # determine the highest category/pheno/evidence grade for this drug_class
                summary_entry.summarise_rules()
                # assign markers with, without rules, and wt markers
                summary_entry.set_markers(flag_core)
                # assign ruleIDs and combo rules
                #TODO: Test combo rule implementation
                # to get the list of possible combo rules to evaluate, we need to extract all 'Combination' rules for this organism
                combo_rules = [r for r in rules if r.get('organism') == summary_entry.organism and r.get('rule type') == 'Combination']
                # then need to further filter to include only combo rules that apply to the drug class we're assessing
                combo_rules = [r for r in combo_rules if summary_entry.drug_class in r.get('drug classes', '')]
                summary_entry.set_ruleIDs_and_combo(combo_rules)
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
                    # determine highest category/pheno/evidence grade for this drug
                    # but take into account the rules for the drug class
                    summary_entry.summarise_rules(master_class_entry)
                    # assign markers
                    summary_entry.set_markers(flag_core)
                    # assign ruleIDs and combo rules
                    combo_rules = [r for r in rules if r.get('organism') == summary_entry.organism and r.get('rule type') == 'Combination']
                    # then need to further filter to include only combo rules that apply to either the drug or class we're assessing
                    combo_rules = [r for r in combo_rules if summary_entry.drug in r.get('drugs', '') or summary_entry.drug_class in r.get('drug classes', '')]
                    summary_entry.set_ruleIDs_and_combo(combo_rules)
                    # add it to our list
                    summary_entry_list.append(summary_entry)
            summary_entry_dict[sample_name] = order_summary_objs(summary_entry_list)
    
    return(summary_entry_dict)