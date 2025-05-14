from amrrules.rules_io import parse_rules_file, download_and_parse_reference_gene_hierarchy
from amrrules.annotator import check_rules, annotate_rule
from amrrules.summariser import prepare_summary, write_output_files
import os, csv

def run(args):
    if args.amr_tool != 'amrfp':
        raise NotImplementedError("Currently only amrfp is supported. Please use amrfp as the AMR tool.")

    if args.amr_tool == 'amrfp':
        # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
        #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
        reference_gene_hierarchy_url = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt"
        amrfp_nodes = download_and_parse_reference_gene_hierarchy(reference_gene_hierarchy_url)

    # let's first deal with the instance of using a single organism across all samples
    rules = parse_rules_file(args.organism, args.amr_tool)

    matched_hits = {}
    unmatched_hits = []
    output_rows = []

    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError("Input file does not contain 'Hierarchy node' column. Please re-run AMRFinderPlus with the --print_node option to ensure this column is in the output file.")
        
        # we need to check if there are multiple samples in this file
        # if there are, we will need to make sure our summary output includes this information
        sample_ids = set()

        row_count = 1
        for row in reader:
            # extract the sample IDs, we will need this later
            # but even if there are multiple sample IDs, it doesn't matter because we're just writing one interpreted output file
            # which is being assessed row by row
            if 'Name' in row:
                sample_ids.add(row.get('Name'))
            matched_rules = check_rules(row, rules, amrfp_nodes)
            if matched_rules is not None:
                # annotate the row with the rule info, based on whether we're using minimal or full annotation
                new_rows = annotate_rule(row, matched_rules, args.annot_opts)
                matched_hits[row_count] = matched_rules
            else:
                new_rows = annotate_rule(row, None, args.annot_opts)
                unmatched_hits.append(row.get('Hierarchy node'))
            row_count += 1
            # add the new rows to the output row list
            output_rows.extend(new_rows)
    
    # okay now we need the unique sample IDs, for the summary output
    sample_ids = list(sample_ids)
    summary_output = prepare_summary(output_rows, rules, sample_ids, args.no_flag_core)

    write_output_files(output_rows, reader, summary_output, args, unmatched_hits, matched_hits)