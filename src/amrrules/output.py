import os
import csv
from amrrules import __version__
from amrrules.utils import required_cols, minimal_columns, full_columns

def write_genotype_report(args, output_rows, unmatched_hits, matched_hits, base_fieldnames):
     # write the output files
    interpreted_output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    #summary_output_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')

    if args.annot_opts == 'minimal':
        interpreted_output_cols = required_cols + minimal_columns
    elif args.annot_opts == 'full':
        interpreted_output_cols = required_cols + minimal_columns + full_columns

    with open(interpreted_output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=base_fieldnames + interpreted_output_cols, delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} hits did not match a rule.")
    print(f"Output written to {interpreted_output_file}.")


def write_genome_report(summary_entry_dict, out_dir, out_prefix):

    # now we want to write out the summary entry report
    # we have all the values we need in each row, under each sample
    summary_output_file = os.path.join(out_dir, out_prefix + '_genome_summary.tsv')
    summary_output_header = ['sample_name', 'drug', 'drug_class', 'category', 'phenotype', 'evidence_grade', 'markers_rule_nonS', 'markers_with_norule', 'markers_S', 'ruleIDs', 'combo_rules', 'organism']
    header_mapping = {
    'sample_name': 'sample',
    'drug': 'drug',
    'drug_class': 'drug class',
    'category': 'clinical category',
    'phenotype': 'phenotype',
    'evidence_grade': 'evidence grade',
    'markers_rule_nonS': 'markers (non-S)',
    'markers_with_norule': 'markers (no rule)',
    'markers_S': 'markers (S)',
    'ruleIDs': 'ruleIDs',
    'combo_rules': 'combo rules',
    'organism': 'organism'   
    }
    csv_header = [header_mapping.get(attr, attr.replace("_", " ").title()) for attr in summary_output_header]

    with open(summary_output_file, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=csv_header, delimiter='\t')
        writer.writeheader()
        for sample, objs in summary_entry_dict.items():
        # Build each row using a dict comprehension, mapping attribute -> CSV header
            rows = [{csv_header[i]: getattr(o, attr, '-') for i, attr in enumerate(summary_output_header)} for o in objs]
            writer.writerows(rows)
    
    print(f"Summary output written to {summary_output_file}.")