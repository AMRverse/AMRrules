from amrrules.rules_io import parse_rules_file, extract_relevant_rules
from amrrules.annotator import check_rules, annotate_rule
from amrrules.summariser import prepare_summary, write_output_files
from amrrules.utils import validate_input_file
from amrrules.resources import ResourceManager as rm
import csv
from importlib import resources

def run(args):

    if args.amr_tool == 'amrfp':
        # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
        #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
        amrfp_nodes = rm().refseq_nodes()

    # We need to validate the input file, and if we've got an organism file, we need to check that all the IDs are in there
    # if only some IDs are missing, in either file, we can still run, but throw a warning to the user
    organism_dict, samples_to_parse = validate_input_file(args.input, args.organism_file, args.organism, args.amr_tool)

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
    output_rows = []

    # now it's time to parse the input file
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')

        row_count = 1
        for row in reader:
            # if we have sample IDs, we know we have the Name column and can get the sample ID
            if samples_to_parse:
                sample_id = row.get('Name')
                # Skip this row if it's not in our sample_to_parse list, as this means we don't have an organism for it
                if sample_id not in samples_to_parse:
                    continue
                # extract the correct organism for this sample ID
                sample_organism = organism_dict.get(sample_id)
            else:
                sample_organism = args.organism
            # extract the relevant rules for this ID, based on its organism
            relevant_rules = extract_relevant_rules(rules, sample_organism)
            matched_rules = check_rules(row, relevant_rules, amrfp_nodes)
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
    summary_output = prepare_summary(output_rows, rules, samples_to_parse, args.no_flag_core)
    #summary_output = None

    write_output_files(output_rows, reader, summary_output, args, unmatched_hits, matched_hits)


def download_resources():
    """
    Download and cache AMRFinderPlus resource files using ResourceManager.
    """
    rm().setup_all_resources()
    print("Resource download complete.")