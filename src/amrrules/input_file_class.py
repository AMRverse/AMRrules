from typing import Any, Dict, Optional
import re

from amrrules.utils import aa_conversion


class InputRow:

    def __init__(self, raw, tool, organism):
        self.raw_row = raw
        self.tool = tool

        # standard fields regardless of tool type
        self.sample_name: Optional[str] = None
        self.gene_symbol: Optional[str] = None
        self.mutation: Optional[str] = None # this will be the formatted AMRrules compliant mutation
        self.amr_type: Optional[str] = None  # type of AMR variant (Gene presence, protein variant, nucl variant etc)
        self.drug_class: Optional[str] = None
        self.drug: Optional[str] = None

        # option to process this row or just skip (eg virulence rows from AMRFP output)
        self.to_process: bool = False

        # amrfp relevant fields(filled by parser)
        self.sample: Optional[str] = None
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
            self.mutation, self.amr_type = self._parse_mutation()
        else:
            self.amr_type = "Gene presence detected"

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

    # Utility / compatibility
    def to_dict(self) -> Dict[str, Any]:
        """Return a dict representation (standardised fields + raw_row under 'raw')."""
        out = {
            "sample": self.sample,
            "gene": self.gene,
            "element_type": self.element_type,
            "element_subtype": self.element_subtype,
            "method": self.method,
            "hierarchy_node": self.hierarchy_node,
            "sequence_accession": self.sequence_accession,
            "hmm_id": self.hmm_id,
            "raw_mutation": self.raw_mutation,
            "mutation_three_letter": self.format_mutation("three_letter"),
            "mutation_one_letter": self.format_mutation("one_letter"),
            "input_type": self.input_type,
        }
        out.update({"raw": dict(self.raw_row)})
        return out

