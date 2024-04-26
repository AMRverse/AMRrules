#!/usr/bin/env python

from argparse import ArgumentParser
import re, os
import csv
from itertools import chain

class GeneRule(object):

    def __init__(self, species, allele, context, drug, expected_pheno):
        self.species = species
        self.allele = allele
        self.context = context
        self.drug = drug
        self.expected_pheno = expected_pheno

class AMRFinderResult(object):

    def __init__(self, allele_id, amr_class, amr_subclass, drug, full_row, name=None):
        self.allele_id = allele_id
        self.amr_class = amr_class
        self.amr_subclass = amr_subclass
        self.drug = str(drug)
        self.name = name
        self.full_row = full_row
        #self.element_type = element_type
        #self.element_subtype = element_subtype
        #self.scope = scope
        #self.sequence_name = sequence_name
        #self.protein_onwards = protein_onwards
        #self.method_onwards = method_onwards

        # Assigned during rule processing
        self.expected_pheno = str()
        self.context = str()
    

class OrganismAwareReport(object):
    
    pass

def get_arguments():
    parser = ArgumentParser(description='Parse AMRFinderPlus files with organism specific rules.')
    
    parser.add_argument('--reports', nargs='+', type=str, required=True, help='One or more AMRFinderPlus results files (should all belong to the same species).')
    parser.add_argument('--species', required=False, default='s__Klebsiella pneumoniae', type=str, help='Species of the genomes in the input files (default is s__Klebsiella pneumoniae)')
    parser.add_argument('--organism_rules', required=True, type=str, help='Organism-specific rule set table, used to generate the annotated AMR report (all genomes included in a single report).')
    parser.add_argument('--drug_dictionary', required=False, default='./example_dict_kleb/Kleb_local_dict.tsv', help='Path to drug dictionary that matches allele names with drug classes. Current default is the temporary dictionary in this repo, "./example_dict_kleb/Kleb_local_dict.tsv')
    parser.add_argument('--output', required=True, type=str, help='Name for output file.')

    return parser.parse_args()

def create_drug_list(drug_gene_file):

    # initialise dictionary, key=allele, value=drug name
    drug_dict = {}

    with open(drug_gene_file, 'r') as drug_genes:
        drug_genes_reader = csv.DictReader(drug_genes, delimiter='\t')
        for row in drug_genes_reader:
            if row["#Allele"] == "NA":
                drug_dict[row["Gene family"]] = str.lower(row["Drug"])
            else:
                drug_dict[row["#Allele"]] = str.lower(row["Drug"])

    return drug_dict

def create_rule_list(rule_infile):
    rule_list = []
    with open(rule_infile, 'r') as rule_file:
        rule_file_reader = csv.DictReader(rule_file, delimiter='\t')
        for row in rule_file_reader:
            new_rule = GeneRule(row["organism"], row["gene"], row["context"], row["drug"], row["category"])
            rule_list.append(new_rule)
    return rule_list

def parse_amr_report(report_file, drug_dict):

    amrfinder_report_lines = []
    
    with open(report_file, 'r') as report:
        report_reader = csv.DictReader(report, delimiter='\t')
        for row in report_reader:
            # only parse the line if the element_type is AMR
            if row["Element type"] == "AMR":
                node_allele = row["Hierarchy node"]
                try:
                    gene_drug = drug_dict[node_allele]
                except KeyError:
                    gene_drug = ''
                try:
                    name_id = row["Name"]
                except KeyError:
                    name_id = None
                amrfinder_result = AMRFinderResult(node_allele, row["Class"], row["Subclass"],
                                                   gene_drug, row, name=name_id)
                amrfinder_report_lines.append(amrfinder_result)

    return amrfinder_report_lines

def determine_rules(amrfinder_report_lines, rule_list, sampleID, species):

    # this is our list of result line classes that have the info we want
    output_lines = []

    # extract the rules relevant to the species
    relevant_rules = [rule for rule in rule_list if rule.species == species]

    for amrfinder_result in amrfinder_report_lines:
        # grab the allele allele AMRFinder found
        amrfinder_allele = amrfinder_result.allele_id
        # okay so how do we do this? Do we take the allele and check that it's in our allele rule list?
        # how are we dealing with cases where we have like blaSHV*? because we want to include all alleles, which assumes then some kind of regex matching or substring in x?
        # allele will be blaSHV-28, we want to match blaSHV*
        for rule in relevant_rules:
            # see if there's a matching rule for that allele
            search_value = re.search(rule.allele, amrfinder_allele)
            # this will return a value if there's something, otherwise it will be None and this won't evaluate
            if search_value:
                # so now we've found something, what do we want to extract? we want to extract the expected phenotype to add to the table, for that rule
                #expanded_result = amrfinder_result.add_to_result(rule.expected_pheno, rule.drug, sampleID, rule.context)

                amrfinder_result.expected_pheno = rule.expected_pheno
                amrfinder_result.drug = rule.drug
                amrfinder_result.name = sampleID
                amrfinder_result.context = rule.context

                # now escape this for loop!!
                break
            # if we're at the final rule, and still no search result then add an empty version, as there is no rule for this allele call
            if (rule_list.index(rule) + 1) == len(rule_list) and not search_value:
                amrfinder_result.name = sampleID
        # add it to our new list
        output_lines.append(amrfinder_result)

    return output_lines

def write_output(output_lines, out_file, species, version):

    with open(out_file, "w") as out:
        # add name to the header if it exists, otherwise don't bother
        header = ['Name', 'Protein identifier', 'Contig id', 'Start', 'Stop', 'Strand', 'Gene symbol', 'Sequence name', 
                  'Species interpretation', 'Context', 'Org interpretation', 'Drug', 'Scope', 'Element type', 'Element subtype', 'Class', 'Subclass',
                  'Method', 'Target length', 'Reference sequence length', '% Coverage of reference sequence', '% Identity to reference sequence',
                  'Alignment length', 'Accession of closest sequence', 'Name of closest sequence', 'HMM id', 'HMM description', 'Hierarchy node']
        out.write('\t'.join(header) + '\n')

        # get the value for species and version for the 'Species interpretation' column
        species_interp = species + '; ' + version

        for out_line in output_lines:
            # correctly format the wt resistant/susceptible codes to match poster
            if out_line.expected_pheno == 'wt resistant':
                expected_pheno = 'wt (R)'
            elif out_line.expected_pheno == 'wt susceptible':
                expected_pheno = 'wt (S)'
            else:
                expected_pheno = out_line.expected_pheno
            
            amrfinder_first_sect = [out_line.full_row['Protein identifier'], out_line.full_row['Contig id'], out_line.full_row['Start'], out_line.full_row['Stop'], out_line.full_row['Strand'], out_line.full_row['Gene symbol'], out_line.full_row['Sequence name']]
            
            org_inter_sect = [species_interp, out_line.context, expected_pheno, out_line.drug]
            
            amrfinder_remaining_sect = [out_line.full_row['Scope'], out_line.full_row['Element type'], out_line.full_row['Element subtype'], out_line.full_row['Class'], out_line.full_row['Subclass'], out_line.full_row['Method'], out_line.full_row['Target length'], out_line.full_row['Reference sequence length'], out_line.full_row['% Coverage of reference sequence'], out_line.full_row['% Identity to reference sequence'], out_line.full_row['Alignment length'], out_line.full_row['Accession of closest sequence'], out_line.full_row['Name of closest sequence'], out_line.full_row['HMM id'], out_line.full_row['HMM description'], out_line.full_row['Hierarchy node']]
            
            final_line = [out_line.name] + amrfinder_first_sect + org_inter_sect + amrfinder_remaining_sect

            out.write('\t'.join(final_line) + '\n')

def organism_aware_report(output_lines, local_drug_list):

    # get the list of drugs we care about, match to their AMRFinder name if needed
    #key: drug name to write in report, value: AMR drug name, if needed, otherwise empty string
    drugs_to_report = []
    drug_amrfinder_name_conversion = {}
    with open(local_drug_list, 'r') as in_file:
        header = 0
        for line in in_file:
            if header == 0:
                header += 1
            else:
                fields = line.strip().split('\t')
                drugs_to_report.append(fields[0])
                if len(fields) == 2:
                    drug_amrfinder_name_conversion[fields[1]] = fields[0]
    
    for entry in output_lines:
        if entry.drug == '' or entry.drug == 'multiple drug':
            drug_interest = entry.amr_subclass
            # check to see the conversion 
        else:
            drug_interest = entry.drug

    return drug_interest

def main():

    args = get_arguments()

    # parse the drug dictionary
    drug_dict = create_drug_list(args.drug_dictionary)

    # parse the organism rule file a list. Each item in the list is an object with the relevant details
    rule_list = create_rule_list(args.organism_rules)


    # to get to the wt(R) thing we want, we're going to need to, for each gene that has a WT R, make a note that this
    #drug is an expected pheno. And then when we get to a call that doesn't have a rule attached to it, check if it's a drug with a known wt(R) expected pheno? but this will depend on the order things happen in
    # kind of need to draw out the class hierachy I think to make this make more sense

    full_output_entries = []

    for report in args.reports:
        amrfinder_results = parse_amr_report(report, drug_dict)
        # if we have a name value in the amrfinder report, use that
        # otherwise just use the name of the report file in the 'Sample' column
        if amrfinder_results[0].name:
            sampleID = amrfinder_results[0].name
        else:
            sampleID = os.path.basename(report)
        output_entries = determine_rules(amrfinder_results, rule_list, sampleID, args.species)

        # sort the lines
        sorted_entries = dict()
        for entry in output_entries:
            # group calls by drug so we can update the expected phenotype as needed
            if entry.drug not in sorted_entries:
                sorted_entries[entry.drug] = list()
            sorted_entries[entry.drug].append(entry)

        # TODO(JH): check out enums
        # check each drug and update to be wt resistant if there is a core gene present
        for drug_name, entries in sorted_entries.items():
            has_wt_r = False
            # Other rules?
            for entry in entries:
                has_wt_r |= entry.expected_pheno == 'wt resistant'

            # Second pass, update information in other entries as required
            for entry in entries:
                if has_wt_r:
                    entry.expected_pheno = 'wt resistant'

        full_output_entries = full_output_entries + output_entries
        
    # now write out the output into a single file
    #TODO: encode version number correctly 
    write_output(full_output_entries, args.output, args.species, 'v1.1')

if __name__ == '__main__':
    main()
