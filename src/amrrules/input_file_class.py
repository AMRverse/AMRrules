from typing import Any, Dict, Optional
import re
from amrrules import __version__
from amrrules.utils import aa_conversion, minimal_columns, full_columns


class Genotype:

    def __init__(self, raw, tool, organism_dict):
        self.raw_row = raw
        self.annotated_rows: Optional[Any] = None # will populate with the anntoated version of the row after rule matching
        self.tool = tool

        # standard fields regardless of tool type
        self.sample_name: Optional[str] = None
        self.gene_symbol: Optional[str] = None
        self.mutation: Optional[str] = None # this will be the formatted AMRrules compliant mutation
        self.variation_type: Optional[str] = None  # type of AMR variant (Gene presence, protein variant, nucl variant etc)
        self.matched_rules: Optional[Any] = None  # will be filled with matched rules

        # option to process this row or just skip (eg virulence rows from AMRFP output)
        self.to_process: bool = False

        # amrfp relevant fields(filled by parser)
        self.nodeID: Optional[str] = None
        self.method: Optional[str] = None
        self.closest_acc: Optional[str] = None
        self.hmm_acc: Optional[str] = None
        self.amrfp_class: Optional[str] = None
        self.amrfp_subclass: Optional[str] = None

        # parse on construction
        if self.tool == "amrfp":
            self._parse_amrfp()
        
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
    def _parse_amrfp(self):
        r = self.raw_row
        element_type = r.get("Element type") or r.get("Type")
        # only process AMR rows
        if element_type != "AMR":
            self.to_process = False
            return
        else:
            self.to_process = True
        
        # get the sample name, but only if the column exists
        if "Name" in r.keys():
            self.sample_name = r.get("Name")
        self.gene_symbol = (r.get("Gene symbol") or r.get("Element symbol"))
        self.method = r.get("Method")
        self.nodeID = r.get("Hierarchy node")
        self.closest_acc = r.get("Accession of closest sequence") or r.get("Closest reference accession")
        self.hmm_acc = r.get("HMM id")

        # if our method is pointX or pointP, we need to extract the actual mutation
        # and convert to AMRrules syntax
        if self.method in ["POINTX", "POINTP", "POINTN"]:
            self.mutation, self.variation_type = self._parse_mutation()
        else:
            self.variation_type = "Gene presence detected"

    def _parse_mutation(self):

        gene_symbol, mutation = self.gene_symbol.rsplit("_", 1)

        # this means it is a protein mutation
        if self.method in ["POINTX", "POINTP"]:
            # extract the relevant parts of the mutation
            pattern = re.compile(r"(\D+)(\d+)(\D+)")
            ref, pos, alt = pattern.match(mutation).groups()
            # convert the single letter AA code to the 3 letter code
            # note that we need to determine if we've got a simple substitution of ref to alt
            # or do we have a deletion or an insertion?
            # gyrA_S83L -> p.Ser83Leu, this is a substitution
            # penA_D346DD -> p.345_346insAsp
            # okay so if there are two characters in alt, then we have an insertion
            if len(alt) > 1 and alt != 'STOP':
                # then we have an insertion, and the inserted AA is the second character of alt
                alt = aa_conversion.get(alt[1])
                # our coordinates are the original pos, and original - 1
                pos_coords = str(int(pos) - 1) + "_" + str(pos)
                return(f"p.{pos_coords}ins{alt}", "Protein variant detected")

            else:
                ref = aa_conversion.get(ref)
                alt = aa_conversion.get(alt)
                return(f"p.{ref}{pos}{alt}", "Protein variant detected")
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

    def _get_final_matches(self, matching_rules):
        if self.variation_type == 'Gene presence detected':
            return matching_rules
        elif self.variation_type in ['Protein variant detected', 'Nucleotide variant detected', 'Promoter variant detected']:
            final_matching_rules = []
            # now we need to check the mutation, extracting any matching rules
            for rule in matching_rules:
                if rule['mutation'] == self.mutation:
                    final_matching_rules.append(rule)
            return final_matching_rules

    def find_matching_rules(self, rules, amrfp_nodes):

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

        #TODO: HMM accession check

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


# New helper class: builds a summary-compatible dict from a ResultRow
class SummaryRow:
    def __init__(self, result_row: ResultRow, card_conversion: Dict = None, card_drug_map: Dict = None):
        """
        result_row: ResultRow instance (annotated_row expected to be set)
        card_conversion: mapping used to convert AMRFP Subclass -> CARD drug/class (optional)
        card_drug_map: mapping from CARD drug -> drug class (optional)
        """
        self.result_row = result_row
        self.card_conversion = card_conversion or {}
        self.card_drug_map = card_drug_map or {}

    def _assign_drug_info(self, card_amrfp_conversion: Dict = None, card_drug_map: Dict = None):
        """
        Assign drug and drug class information
        """
        # if we have a matched rule, extract that info
        # note that if we get drug, we won't have a class assigned so need to fill that in
        if self.matched_rule:
            self.drug = self.matched_rule.get('drug', '-')
            if self.drug != '-':
                # get the drug class from card
                self.drug_class = card_drug_map.get(self.drug, '-')
            else:
                self.drug_class = self.matched_rule.get('drug class', '-')
        # if we DON'T have a matched rule, we need to assign it from the AMRFP information
        # using our conversion dictionary
        else:
            subclasses = self.genotype.amrfp_subclass.split('/') if '/' in self.genotype.amrfp_subclass else [self.genotype.amrfp_subclass]
            # okay if subclasses has more than one entry.... omg we need to duplicate this for each subclass as we'll have different drugs
            # so this is ONLY for instances where we have no matched rule, and it's an AMR marker
            # let's save this as a list and deal with it in the summary row section
            self.amr_subclass_drugs = []
            for subclass in subclasses:
                drug = self.card_conversion.get(subclass, {}).get('drug', '-')
                drug_class = self.card_drug_map.get(drug, '-')
                self.amr_subclass_drugs.append((subclass, drug, drug_class))
    
    def to_summary_dict(self):
        """
        Return a dict compatible with the summariser expectations.
        If the ruleID is '-' (no rule) and a Subclass exists, try to infer drug/drug class from card_conversion.
        """
        # use annotated row if available, otherwise fallback to base_row
        row = dict(self.result_row.annotated_row) if self.result_row.annotated_row else dict(self.result_row.base_row)

        # Ensure keys used in summariser exist
        row.setdefault('ruleID', row.get('ruleID', '-'))
        row.setdefault('Subclass', row.get('Subclass', 'NA'))
        row.setdefault('drug', row.get('drug', '-'))
        row.setdefault('drug class', row.get('drug class', '-'))
        row.setdefault('phenotype', row.get('phenotype', '-'))
        row.setdefault('clinical category', row.get('clinical category', '-'))
        row.setdefault('evidence grade', row.get('evidence grade', '-'))
        row.setdefault('organism', row.get('organism', self.result_row.organism or '-'))
        row.setdefault('Name', row.get('Name', self.result_row.sample_name))
        row.setdefault('Gene symbol', row.get('Gene symbol', row.get('Element symbol')))
        row.setdefault('Element symbol', row.get('Element symbol', row.get('Gene symbol')))
        row.setdefault('context', row.get('context', '-'))

        # If no rule assigned and Subclass present, attempt to resolve CARD drug/class
        if (row.get('ruleID') == '-' or row.get('ruleID') is None) and row.get('Subclass') and row.get('Subclass') != 'NA':
            subclasses = row.get('Subclass')
            # handle slash-separated subclasses
            subclasses_list = subclasses.split('/') if '/' in subclasses else [subclasses]
            # pick first resolvable mapping; leave defaults if none found
            for sc in subclasses_list:
                conv = self.card_conversion.get(sc, {})
                resolved_drug = conv.get('drug')
                resolved_class = conv.get('class')
                if resolved_drug and resolved_drug not in ('NA','-'):
                    row['drug'] = resolved_drug
                    # try to fill drug class using provided map if not present
                    row['drug class'] = self.card_drug_map.get(resolved_drug, row.get('drug class', '-'))
                    break
                if (not resolved_drug or resolved_drug in ('NA','-')) and resolved_class and resolved_class not in ('NA','-'):
                    row['drug class'] = resolved_class
                    break

            # If still NA, mark as other markers and set some fields to '-'
            if row['drug'] in ('NA','-') and row['drug class'] in ('NA','-'):
                row['drug'] = 'other markers'
                row['drug class'] = 'other markers'
                row['phenotype'] = row.get('phenotype', '-')
                row['clinical category'] = row.get('clinical category', '-')
                row['evidence grade'] = row.get('evidence grade', '-')

        return row


