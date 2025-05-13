from amrrules.rules_io import parse_rules_file, download_and_parse_reference_gene_hierarchy
from amrrules.annotator import check_rules, annotate_rule
from amrrules.summariser import write_summary
import os, csv

def run(args):
    if args.amr_tool != 'amrfp':
        raise NotImplementedError("Currently only amrfp is supported.")

    nodes = download_and_parse_reference_gene_hierarchy()
    rules = parse_rules_file(args.rules)

    matched_hits = {}
    unmatched_hits = []
    output_rows = []

    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        if 'Hierarchy node' not in reader.fieldnames:
            raise ValueError("Missing 'Hierarchy node' column.")

        for i, row in enumerate(reader, 1):
            matched_rules = check_rules(row, rules, nodes)
            if matched_rules:
                new_rows = annotate_rule(row, matched_rules, args.annot_opts)
                matched_hits[i] = matched_rules
            else:
                new_rows = annotate_rule(row, None, args.annot_opts)
                unmatched_hits.append(row.get('Hierarchy node'))

            output_rows.extend(new_rows)

    _write_outputs(reader.fieldnames, output_rows, args, rules, matched_hits, unmatched_hits)

def _write_outputs(fieldnames, output_rows, args, rules, matched_hits, unmatched_hits):
    output_file = os.path.join(args.output_dir, args.output_prefix + '_interpreted.tsv')
    summary_file = os.path.join(args.output_dir, args.output_prefix + '_summary.tsv')

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames + ['ruleID', 'context', 'drug', 'drug class', 'phenotype', 'clinical category', 'evidence grade'], delimiter='\t')
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"{len(matched_hits)} hits matched a rule and {len(unmatched_hits)} didn't.")
    print(f"Output written to {output_file}")

    summary_rows = write_summary(output_rows, rules)
    with open(summary_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['drug', 'drug class', 'category', 'phenotype', 'evidence grade', 'markers', 'ruleIDs', 'combo rules'], delimiter='\t')
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Summary written to {summary_file}")