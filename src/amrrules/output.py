import os
import csv

def write_genotype_report(args, output_rows, reader, unmatched_hits, matched_hits):
     # write the output files
    interpreted_output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    #summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')

    minimal_columns = ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade', 'version', 'organism']
    full_columns = ['breakpoint', 'breakpoint standard', 'breakpoint condition', 'evidence code', 'evidence limitations', 'PMID', 'rule curation note']
    if args.annot_opts == 'minimal':
        interpreted_output_cols = minimal_columns
    elif args.annot_opts == 'full':
        interpreted_output_cols = minimal_columns + full_columns

    with open(interpreted_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=reader.fieldnames + interpreted_output_cols, delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} hits did not match a rule.")
    print(f"Output written to {interpreted_output_file}.")

def write_genome_report(args, summary_output):

    summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_geome_summary.tsv')
    with open(summary_output_file, 'w', newline='') as f:
        if 'Name' in summary_output[0].keys():
            col_names = ['Name', 'drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers wt', 'markers nwt', 'ruleIDs', 'combo rules', 'organism']
        else:
            col_names = ['drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers wt', 'markers nwt', 'ruleIDs', 'combo rules', 'organism']
        writer = csv.DictWriter(f, fieldnames=col_names, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_output)
    print(f"Summary output written to {summary_output_file}.")