from amrrules.rules_io import parse_rules_file, download_and_parse_reference_gene_hierarchy
from amrrules.annotator import check_rules, annotate_rule
from amrrules.summariser import write_summary
import os, csv

def run(args):
    if args.amr_tool != 'amrfp':
        raise NotImplementedError("Currently only amrfp is supported. Please use amrfp as the AMR tool.")

    if args.amr_tool == 'amrfp':
        # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
        #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
        reference_gene_hierarchy_url = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt"
        amrfp_nodes = download_and_parse_reference_gene_hierarchy(reference_gene_hierarchy_url)

    rules = parse_rules_file(args.rules, args.amr_tool)

    matched_hits = {}
    unmatched_hits = []
    output_rows = []

    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError("Input file does not contain 'Hierarchy node' column. Please re-run AMRFinderPlus with the --print_node option to ensure this column is in the output file.")

        row_count = 1
        for row in reader:
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

    # write the output files
    interpreted_output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')
    with open(interpreted_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames + ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade', 'version'], delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} hits did not match a rule.")
    print(f"Output written to {interpreted_output_file}.")

    summary_output = write_summary(output_rows, rules)

    with open(summary_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers', 'ruleIDs', 'combo rules'], delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_output)
    print(f"Summary output written to {summary_output_file}.")