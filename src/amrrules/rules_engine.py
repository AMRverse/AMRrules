from amrrules.rules_io import parse_rules_file, extract_relevant_rules
from amrrules.summariser import create_summary_dict
from amrrules.utils import check_sample_ids, validate_amrfp_file, get_organisms
from amrrules.output import write_genotype_report, write_genome_report
from amrrules.resources import ResourceManager as rm
from amrrules.genotype_parser import GenoResult, Genotype
import csv
from importlib import resources
from collections import defaultdict

def run(args):

    # extract all the rules relevant to the organisms we're processing
    if args.organism_file:
        organism_dict, skipped_samples = get_organisms(args.organism_file)
        samples_with_org = set(organism_dict.keys())
    else:
        organism_dict = {'': args.organism}
        samples_with_org = None
        skipped_samples = None
    
    if args.amr_tool == 'amrfp':
        # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
        #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
        amrfp_nodes = rm().refseq_nodes()
        # check the input file has the Hierarchy node column, and if an organism file is included, that the first column is Name
        samples_to_parse = validate_amrfp_file(args.input, multi_entry=bool(args.organism_file))
        
        # get the AMRFP to CARD conversion mapping for later use - we only want to do this once
        # so do it here and pass this to where it's needed
        card_amrfp_conversion = rm().get_amrfp_card_conversion() # get the AMRFP to CARD conversion mapping
        card_drug_map = rm().get_card_drug_class_map() # get CARD drugs and their associated classes
    
    # if the is a multi-entry file, we need to check that all our sampleIDs are in the organism file
    # will raise an error if any are missing
    # will raise a warning if there are samples in the org file but aren't in the input file
    if args.organism_file:
        check_sample_ids(samples_with_org, samples_to_parse, skipped_samples)

    # collate the required rules for the organisms we need to parse
    rule_files = []
    # open the rules key file and get the organism name
    key_file_path = resources.files("amrrules.rules").joinpath("rule_key_file.tsv")
    with open(key_file_path, 'r') as key_file:
        for row in key_file:
            # split the row into the organism and rules file
            organism, rules_filename = row.strip().split('\t')
            # if it's an organism we're interested in, add it to the list of rule files to parse
            if organism in set(organism_dict.values()):
                rule_files.append(rules_filename)

    # parse the rule files
    rules = parse_rules_file(rule_files)

    matched_hits = {}
    unmatched_hits = []
    genotype_rows = []
    # now it's time to parse the input file, which we have validated to check that it has
    # the columns we need. Each row will be parsed into an InputRow object
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        base_fieldnames = reader.fieldnames.copy()
        row_count = 1
        for row in reader:
            if args.sample_id:
                row_to_process = GenoResult(row, args.amr_tool, organism_dict, sample_name=args.sample_id)
            else:
                row_to_process = GenoResult(row, args.amr_tool, organism_dict)
            # if this row belongs to a sample we should skip, update the to_process attribute to False
            if skipped_samples and row_to_process.sample_name in skipped_samples:
                row_to_process.to_process = False
            # we only want to find matched rules for a row if it's relevant for AMR, so check this value first
            # also make sure it's not a row belonging to a sample we should skip
            if row_to_process.to_process:                
                # extract the relevant rules for this ID, based on its organism
                relevant_rules = extract_relevant_rules(rules, row_to_process.organism)
                # determine if there's a matching rule for this row (this sets row_to_process.matched_rules)
                row_to_process.find_matching_rules(relevant_rules, amrfp_nodes)
            
            row_to_process.annotate_row(args.annot_opts)

            # track matched / unmatched hits for reporting
            # create a result row for each matched rule, as we need to duplicate rows in output
            # if they have multiple matching rules
            if row_to_process.matched_rules:
                matched_hits[row_count] = row_to_process.matched_rules
            # if there's no matching rule, or it's a row we don't process, still create a ResultRow, but it has no rule
            else:
                unmatched_hits.append(row.get('Hierarchy node'))

            # keep Genotype objects in case we need them later
            genotype_rows.append(row_to_process)
            row_count += 1
    
    # get all the output rows together into a single list
    genotype_output_rows = []
    for g in genotype_rows:
        genotype_output_rows.extend(g.annotated_row)

    # now write out the interpreted genotype report, which annotates each row with the rule info
    write_genotype_report(args, genotype_output_rows, unmatched_hits, matched_hits, base_fieldnames)

    # we now want to create one object per rule/AMRFP subclass, so that we can summarise by drug or drug class.
    genotype_objects = []
    for g in genotype_rows:
        if g.to_process:
            if g.matched_rules:
                for rule in g.matched_rules:
                    geno_obj = Genotype.from_result_row(g, card_map=card_drug_map, rule=rule)
                    genotype_objects.append(geno_obj)
            else:
                # extract the subclasses and split as needed
                g_subclasses = g.amrfp_subclass.split('/')
                for subclass in g_subclasses:
                    geno_obj = Genotype.from_result_row(g, card_amrfp=card_amrfp_conversion, amrfp_subclass=subclass, no_rule_interp=args.no_rule_interpretation)
                    genotype_objects.append(geno_obj)

    # now we want to group all of these objects by sample ID (if we have multiple samples)
    # because we need to summarise per genome
    # then we want to group by drug, or drug class if drug is '-', in each sample
    grouped_by_sample = defaultdict(list)
    for geno_obj in genotype_objects:
        grouped_by_sample[geno_obj.sample_name].append(geno_obj)


    summary_entry_dict = create_summary_dict(grouped_by_sample, rules)
    
    write_genome_report(summary_entry_dict, args.output_dir, args.output_prefix)


def download_resources():
    """
    Download and cache the AMRFP and CARD database files required.
    """
    rm().setup_all_resources()
    print("Resource download complete.")