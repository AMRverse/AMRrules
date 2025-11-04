import os
import csv
from amrrules import __version__
from amrrules.utils import minimal_columns, full_columns

def annotate_row(row, rules, annot_opts, version=__version__):
    if rules is None:
        # if we didn't find a matching rule, then we need to add new columns for each of the options but using '-' as the value
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = '-'
        elif annot_opts == 'full':
            for col in minimal_columns + full_columns:
                row[col] = '-'
        row['version'] = version
        row['organism'] = '-'
        return [row]
    # if we found multiple rules, we want to return each rule as its own row in the output file
    if len(rules) > 1:
        # we need to create a new row for each rule
        output_rows = []
        for rule in rules:
            new_row = row.copy()
            if annot_opts == 'minimal':
                for col in minimal_columns:
                    new_row[col] = rule.get(col)
            elif annot_opts == 'full':
                for col in minimal_columns + full_columns:
                    new_row[col] = rule.get(col)
            new_row['version'] = version
            output_rows.append(new_row)
        return output_rows
    # if we found a single rule, we want to annotate the row with the rule info
    else:
        if annot_opts == 'minimal':
            for col in minimal_columns:
                row[col] = rules[0].get(col)
        elif annot_opts == 'full':
            for col in minimal_columns + full_columns:
                row[col] = rules[0].get(col)
        row['version'] = version
        return [row]

def write_genotype_report(args, output_rows, unmatched_hits, matched_hits, base_fieldnames):
     # write the output files
    interpreted_output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    #summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')

    if args.annot_opts == 'minimal':
        interpreted_output_cols = minimal_columns
    elif args.annot_opts == 'full':
        interpreted_output_cols = minimal_columns + full_columns

    with open(interpreted_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=base_fieldnames + interpreted_output_cols, delimiter='\t')
        writer.writeheader()
        for row in output_rows:
            writer.writerow(row)
    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} hits did not match a rule.")
    print(f"Output written to {interpreted_output_file}.")

def write_genome_report(args, summary_output):

    summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_genome_summary.tsv')
    col_names = ['drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers (with rule)', 'markers (no rule)', 'wt markers', 'ruleIDs', 'combo rules', 'organism']
    with open(summary_output_file, 'w', newline='') as f:
        if 'Name' in summary_output[0].keys():
            col_names = ['Name'] + col_names
        writer = csv.DictWriter(f, fieldnames=col_names, delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_output)
    print(f"Summary output written to {summary_output_file}.")