from typing import Any, Optional
import re
from amrrules import __version__
from amrrules.utils import aa_conversion, minimal_columns, full_columns


class GenoResult:

    def __init__(self, raw, tool, organism_dict, print_non_amr, full_disrupt, sample_name=None):
        self.raw_row = raw
        self.annotated_rows: Optional[Any] = None # will populate with the anntoated version of the row after rule matching
        self.tool = tool
        self.sample_name = sample_name

        # standard fields regardless of tool type
        self.gene_symbol: Optional[str] = None
        self.marker_amrrules: Optional[str] = None  # this will be the formatted AMRrules compliant marker, format gene:mutation (if mutation exists)
        self.mutation: Optional[str] = None # this will be the formatted AMRrules compliant mutation
        self.variation_type: Optional[str] = None  # type of AMR variant (Gene presence, protein variant, nucl variant etc)
        self.matched_rules: Optional[Any] = None  # will be filled with matched rules

        # option to process this row or just skip (eg virulence rows from AMRFP output)
        self.to_process: bool = False
        # option to print this row into the interpreted output, default is not to print
        self.print_row: bool = False

        # amrfp relevant fields(filled by parser)
        self.nodeID: Optional[str] = None
        self.subtype: Optional[str] = None
        self.method: Optional[str] = None
        self.closest_acc: Optional[str] = None
        self.hmm_acc: Optional[str] = None
        self.amrfp_class: Optional[str] = None
        self.amrfp_subclass: Optional[str] = None

        # parse on construction
        if self.tool == "amrfp":
            self._parse_amrfp(print_non_amr, full_disrupt)
        
        #TODO: implement parsing of other input types
        #elif self.input_type == "card":
        #    self._parse_card()
        #elif self.input_type in ("harmonized", "hamronized", "harmonised"):
        #    self._parse_harmonized()
        
        # get the organism info for this sample
        if len(organism_dict.keys()) == 1:
            # then we have a single organism regardless of row
            self.organism = organism_dict.get('')
        else:
            # use the sample name to get the organism
            self.organism = organism_dict.get(self.sample_name)

    # Parsing helpers for specific input types
    def _parse_amrfp(self, print_non_amr, full_disrupt):
        r = self.raw_row
        element_type = r.get("Element type") or r.get("Type")
        # only process AMR rows
        if element_type != "AMR":
            self.to_process = False
            # default is to skip non-AMR rows, but these can be put back in with user optio
            if print_non_amr:
                self.print_row = True
            return
        else:
            self.to_process = True
            self.print_row = True
        
        # get the sample name, but only if the column exists
        # otherwise we will use the sample name provided by the user
        if "Name" in r.keys() and not self.sample_name:
            self.sample_name = r.get("Name")
        # if we weren't provide a sample name, and also we don't have one
        # from the input file, we will use the default name "sample" instead
        if "Name" not in r.keys() and not self.sample_name:
            self.sample_name = "sample"
        self.gene_symbol = (r.get("Gene symbol") or r.get("Element symbol"))
        self.subtype = r.get("Subtype") or r.get("Element subtype")
        self.method = r.get("Method")
        self.nodeID = r.get("Hierarchy node")
        self.closest_acc = r.get("Accession of closest sequence") or r.get("Closest reference accession")
        self.hmm_acc = r.get("HMM id") or r.get("HMM accession")
        self.amrfp_class = r.get("Class")
        self.amrfp_subclass = r.get("Subclass")

        # if our method is pointX or pointP, we need to extract the actual mutation
        # and convert to AMRrules syntax
        if self.method in ["POINTX", "POINTP", "POINTN"]:
            self.mutation, self.variation_type = self._parse_mutation()
        # if there's an internal stop, or the gene is a partial hit, then it's inactivated
        # ditto for POINT_DISRUPT
        elif self.method == "INTERNAL_STOP" or self.method.startswith("PARTIAL"):
            self.variation_type = "Inactivating mutation detected"
            self.mutation = '-' # set mutation to '-' for the summary marker later
        else:
            self.variation_type = "Gene presence detected"
        
        # create the AMRrules compliant marker
        self.marker_amrrules = self._create_amrrules_marker(full_disrupt)

    def _parse_mutation(self):

        gene_symbol, mutation = self.gene_symbol.rsplit("_", 1)
        # set the amrrules marker to the gene symbol for now
        self.marker_amrrules = gene_symbol

        # this means it is a protein mutation
        if self.method in ["POINTX", "POINTP"]:
            # variation type is protein variant, unless subtype is POINT_DISRUPT
            # which means the variation type is inactivation mutation
            if self.subtype == "POINT_DISRUPT":
                variation_type = "Inactivating mutation detected"
            else:
                variation_type = "Protein variant detected"
            # extract the relevant parts of the mutation
            pattern = re.compile(r"(\D+)(\d+)(\D+)")
            ref, pos, alt = pattern.match(mutation).groups()
            # convert the single letter AA code to the 3 letter code
            # note that we need to determine if we've got a simple substitution of ref to alt
            # or do we have a deletion or an insertion, or a frameshift?
            # gyrA_S83L -> p.Ser83Leu, this is a substitution
            # penA_D346DD -> p.345_346insAsp, this is an insertion
            # ompK35_E42RfsTer47 -> p.Glu42ArgfsTer47, this is an inactivating frameshift
            # ompK35_Y36Ter -> p.Tyr36Ter, this is an inactivating frameshift

            # okay so if there are two or more characters in alt, and no Ter, then we have an insertion
            if len(alt) > 1 and alt != 'STOP' and 'Ter' not in alt:
                # then we have an insertion, and the inserted AAs is the second character onwards of alt
                # need to get all the inserted AAs converted
                alt_string = ''
                for char in alt[1:]:
                    alt_string += aa_conversion.get(char)
                # our coordinates are the original pos, and original - 1
                pos_coords = str(int(pos) - 1) + "_" + str(pos)
                return(f"p.{pos_coords}ins{alt_string}", variation_type)
            if len(alt) > 1 and 'Ter' in alt:
                # then we have a frameshift leading to a stop codon
                # if the start of alt is Ter, then it's simply ref pos *
                if alt.startswith('Ter'):
                    return(f"p.{aa_conversion.get(ref)}{pos}Ter", variation_type)
                else:
                    # otherwise we need to get the first aa and convert it
                    alt = aa_conversion.get(alt[0])
                    # then grab the number after Ter
                    fs_pos = re.search(r'Ter(\d+)', mutation).group(1)
                    return(f"p.{aa_conversion.get(ref)}{pos}{alt}fsTer{fs_pos}", variation_type)
            # this is the option if we have a simple conversion of one aa to another
            else:
                ref = aa_conversion.get(ref)
                alt = aa_conversion.get(alt)
                return(f"p.{ref}{pos}{alt}", variation_type)
        
        elif self.method == "POINTN":
            # we need to extract the relevant parts, this will be different because we may have promoter mutations
            pattern = re.compile(r'^([A-Za-z]+)(-?\d+)([A-Za-z]+)$')
            ref, pos, alt = pattern.match(mutation).groups()
            if '-' in pos:
                mutation_type = "Promoter variant detected"
            else:
                mutation_type = "Nucleotide variant detected"
            # if there's a '-' in the position, then this is a promoter mutation
            # if 'del' is in mutation, then we need to convert to c.PosNTdel.
            if alt == 'del':
                return(f"c.{pos}{ref}del", mutation_type)
            # otherwise it's more like 23S_G2032T -> c.2032G>T, with a - if it's in the promoter.
            else:
                return(f"c.{pos}{ref}>{alt}", mutation_type)
    
    def _create_amrrules_marker(self, full_disrupt):
        # if variation type is gene presence, then just return the gene symbol
        if self.variation_type == "Gene presence detected":
            return self.gene_symbol
        # if the variation type is an inactivating mutation, return gene:-, where - indicates inactivation
        #elif self.variation_type == "Inactivating mutation detected":
        #    return f"{self.gene_symbol}:-"
        # otherwise return gene:mutation, which will be '-' if it's an inactivating mutation with no specifics
        else:
            if self.marker_amrrules:
                # if we have a POINT_DISRUPT and the user has set full-disrupt option, then show the full mutation
                if self.subtype == "POINT_DISRUPT" and full_disrupt:
                    return f"{self.marker_amrrules}:{self.mutation}"
                # otherwise return just a '-' to indicate it's inactivated
                elif self.subtype == "POINT_DISRUPT" and not full_disrupt:
                    return f"{self.marker_amrrules}:-"
                # all other cases return gene:mutation
                else:
                    return f"{self.marker_amrrules}:{self.mutation}"
            else:
                return f"{self.gene_symbol}:{self.mutation}"

    def _get_final_matches(self, matching_rules, guideline_pref = None):
        #TODO MAKE THIS ITS OWN FUNCTION
        # we have multiple rules that match, and if there's a guideline preference
        # we need to pick one
        #if guideline_pref:
        #    for rule in matching_rules:
                    # extract the first word inside breakpoint standard (which will be either EUCAST or CLSI)
        #            guideline_for_rule = rule['breakpoint standard'].split()[0]
        if self.variation_type == 'Gene presence detected':
            return matching_rules
        elif self.variation_type in ['Protein variant detected', 'Nucleotide variant detected', 'Promoter variant detected']:
            final_matching_rules = []
            # now we need to check the mutation, extracting any matching rules
            for rule in matching_rules:
                if rule['mutation'] == self.mutation:
                    final_matching_rules.append(rule)
            return final_matching_rules

    def find_matching_rules(self, rules, amrfp_nodes, guideline_pref = None):

        # select rules that match our variation type
        rules_to_check = []
        for rule in rules:
                if rule['variation type'] == self.variation_type:
                    rules_to_check.append(rule)

        # First we're going to check for the nodeID, and if we have one or matches, we we return that
        matching_rules = [rule for rule in rules_to_check if rule.get('nodeID') == self.nodeID]
        if len(matching_rules) > 0:
            self.matched_rules = self._get_final_matches(matching_rules)
            return

        # Okay so nothing matched directly to the nodeID, or we would've returned out of the function. 
        # So now we need to check if there's a parent node that matches our nodeID
        # Note we don't need to check for variation type here, because we're not going to need to go up the hierarchy
        # for that type of rule
        parent_node = amrfp_nodes.get(self.nodeID)
        while parent_node is not None and parent_node != 'AMR':
            matching_rules = [rule for rule in rules_to_check if rule.get('nodeID') == parent_node]
            if len(matching_rules) > 0:
                self.matched_rules = self._get_final_matches(matching_rules)
                return
            parent_node = amrfp_nodes.get(parent_node)

        #Okay so using the nodeID didn't work, so now we need to check the sequence accession
        # start with the nucleotide accessions
        matching_rules = [rule for rule in rules_to_check if rule.get('nucleotide accession') == self.closest_acc]
        if len(matching_rules) > 0:
            self.matched_rules = self._get_final_matches(matching_rules)
            return
        # then check the protein accessions
        matching_rules = [rule for rule in rules_to_check if rule.get('protein accession') == self.closest_acc]
        if len(matching_rules) > 0:
            self.matched_rules = self._get_final_matches(matching_rules)
            return

        #HMM accession check
        matching_rules = [rule for rule in rules_to_check if rule.get('HMM accession') == self.hmm_acc]
        if len(matching_rules) > 0:
            self.matched_rules = self._get_final_matches(matching_rules)
            return

        # if nothing matched, then we return and the value stays the default which is None
        return
    
    def annotate_row(self, annot_opts: str):
        """
        Annotate the base_row using the matched_rule(s) and store in annotated_row.

        Parameters:
            annot_opts (str): Either 'minimal' or 'full'.
                - 'minimal': Only minimal_columns are annotated.
                - 'full': Both minimal_columns and full_columns are annotated.

        Returns:
            List[Dict]: A list of dictionaries containing the annotated row(s).
        """
        base_row = self.raw_row
        annotated_rows = []
        cols = minimal_columns if annot_opts == 'minimal' else minimal_columns + full_columns
        # fill 'variation type', 'gene' and 'mutation' columns based on parsed info, regardless 
        # of whether there's a matched rule or not
        if self.to_process:
            base_row['variation type'] = self.variation_type
            base_row['gene'] = self.marker_amrrules.split(':')[0]  # just the gene part
            base_row['mutation'] = self.mutation if self.mutation else '-'
        else:
            base_row['variation type'] = 'Non-AMR element'
            base_row['gene'] = self.gene_symbol if self.gene_symbol else '-'
            base_row['mutation'] = '-'

        # No matching rules: return single row with '-' for annotation columns
        if not self.matched_rules:
            for col in cols:
                base_row[col] = '-'
            base_row['version'] = __version__
            annotated_rows.append(base_row)
            self.annotated_row = annotated_rows

        # One or more matching rules: create one annotated row per rule
        else:
            if len(self.matched_rules) > 0:
                rules_to_use = self.matched_rules
            else:
                rules_to_use = [self.matched_rule]
            for rule in rules_to_use:
                row = base_row.copy()
                for col in cols:
                    row[col] = rule.get(col, '-')
                row['version'] = __version__
                # prefer organism from rule if present, else from genotype
                row['organism'] = self.organism
                annotated_rows.append(row)

        self.annotated_row = annotated_rows
        return annotated_rows

# we now need to take our genotype objects, and instead group them by drug (or class if no drug specified)
# so each genotype object may have multiple drugs associated with it, regardless of whether it has a matched rule or not

class Genotype(GenoResult):
    def __init__(self, row, tool, organism_dict, matching_rule=None, amrfp_subclass=None):
        super().__init__(row, tool, organism_dict)
        self.matching_rule = matching_rule
        self.amrfp_subclass = amrfp_subclass
        self.drug: Optional[str] = None
        self.drug_class: Optional[str] = None
    
    @classmethod
    def from_result_row(cls, geno_result_obj, card_amrfp=None, card_map=None, rule=None, amrfp_subclass=None, no_rule_interp = None):
        """
        Create a Genotype instance from an existing GenoResult object,
        copying all of its attributes.
        """
        # Create an empty Genotype object without re-parsing anything
        # by bypassing GenoResult.__init__
        new_obj = cls.__new__(cls)

        # Copy all data from the GenoResult instance
        new_obj.__dict__ = dict(geno_result_obj.__dict__)

        # asign the new attributes
        new_obj.rule = rule
        new_obj.amrfp_subclass = amrfp_subclass
        new_obj.has_rule = False

        # Assign drugs if the rule is provided
        # also assign the important values from the rule that are relevant for our summary functions
        if new_obj.rule:
            new_obj._assign_drug_from_rule(card_map)
            new_obj.has_rule = True
            new_obj._assign_rule_attributes(rule)
        elif new_obj.amrfp_subclass:
            new_obj._assign_drug_from_amrfp(card_amrfp)
            new_obj._assign_norule_attributes(no_rule_interp)

        return new_obj

    def _assign_drug_from_rule(self, card_drug_map):
        self.drug = self.rule.get('drug', '-')
        if self.drug != '-':
            # get the drug class from card
            self.drug_class = card_drug_map.get(self.drug, '-')
        else:
            self.drug_class = self.rule.get('drug class', '-')
    
    def _assign_drug_from_amrfp(self, card_amrfp_conversion):
        self.drug = card_amrfp_conversion.get(self.amrfp_subclass).get('drug', '-')
        self.drug_class = card_amrfp_conversion.get(self.amrfp_subclass).get('class', '-')
    
    def _assign_rule_attributes(self, rule):
        # assign other important attributes from the rule for summary purposes
        self.gene_context = rule.get('gene context')
        self.phenotype = rule.get('phenotype')
        self.clinical_category = rule.get('clinical category')
        self.evidence_grade = rule.get('evidence grade')
        self.ruleID = rule.get('ruleID')
    
    def _assign_norule_attributes(self, no_rule_interpretation):
        # assign default values when no rule is matched
        # note that our interpretation depends on user choice
        # AND also if the variation type is "Inactivating mutation detected" but we've got no rule

        # regardless of user choice, the phenotype is always nonwildtype
        self.phenotype = 'nonwildtype'
        # if the variation type is inactivating, and we've got no rule
        # then we treat it as though the gene isn't functional, so therefore is S
        if self.variation_type == "Inactivating mutation detected":
            # the ONLY exception to this is if it's a POINT_DISRUPT, in which case the default interpretation should be R
            # (if nwtR), as these are interruptions in core genes that are likely to cause resistance
            if self.subtype == "POINT_DISRUPT" and no_rule_interpretation == 'nwtR':
                self.clinical_category = 'R'
            else:
                self.clinical_category = 'S'
        # otherwise the gene is considered functional, so we use what the user has selected
        if self.variation_type != "Inactivating mutation detected" and no_rule_interpretation == 'nwtR':
            self.clinical_category = 'R'
        elif no_rule_interpretation == 'nwtS':
            self.clinical_category = 'S'
        
        # regardless evidence is very low and there is no ruleID to set
        self.evidence_grade = 'very low'
        self.ruleID = None