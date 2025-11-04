from typing import Any, Dict, Optional
import re
from amrrules.utils import aa_conversion


class InputRow:

    def __init__(self, raw, tool, organism_dict):
        self.raw_row = raw
        self.tool = tool

        # standard fields regardless of tool type
        self.sample_name: Optional[str] = None
        self.gene_symbol: Optional[str] = None
        self.mutation: Optional[str] = None # this will be the formatted AMRrules compliant mutation
        self.variation_type: Optional[str] = None  # type of AMR variant (Gene presence, protein variant, nucl variant etc)
        self.drug_class: Optional[str] = None
        self.drug: Optional[str] = None
        self.matched_rules: Optional[Any] = None  # will be filled with matched rules

        # option to process this row or just skip (eg virulence rows from AMRFP output)
        self.to_process: bool = False

        # amrfp relevant fields(filled by parser)
        self.nodeID: Optional[str] = None
        self.method: Optional[str] = None
        self.closest_acc: Optional[str] = None
        self.hmm_acc: Optional[str] = None

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


