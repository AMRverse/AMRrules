from amrrules.rules_io import parse_rules_file, extract_relevant_rules
from amrrules.summariser import prepare_summary
from amrrules.utils import check_sample_ids, validate_amrfp_file, get_organisms
from amrrules.output import write_genotype_report, write_genome_report
from amrrules.resources import ResourceManager as rm
from amrrules.input_file_class import Genotype, ResultRow, SummaryRow
import csv
from importlib import resources

def run(args):

    # extract all the rules relevant to the organisms we're processing
    if args.organism_file:
        organism_dict = get_organisms(args.organism_file)
        samples_with_org = set(organism_dict.keys())
    else:
        organism_dict = {'': args.organism}
        samples_with_org = None
    
    if args.amr_tool == 'amrfp':
        # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
        #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
        amrfp_nodes = rm().refseq_nodes()
        # check the input file has the Hierarchy node column, and if an organism file is included, that the first column is Name
        samples_to_parse = validate_amrfp_file(args.input, multi_entry=bool(args.organism_file))
    
    # if the is a multi-entry file, we need to check that all our sampleIDs are in the organism file
    # will raise an error if any are missing
    # will raise a warning if there are samples in the org file but aren't in the input file
    if args.organism_file:
        check_sample_ids(samples_with_org, samples_to_parse)

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
    # TODO: consider convert rules to Rule objects? Not sure how that helps us though
    #rules = [Rule(r) for r in rules]

    # We need to validate the input file, and if we've got an organism file, we need to check that all the IDs are in there
    # if only some IDs are missing, in either file, we can still run, but throw a warning to the user

    matched_hits = {}
    unmatched_hits = []
    output_rows = []

    genotype_rows = []
    result_row_objs = []
    # collect rows prepared for summariser (dicts)
    summary_input_rows = []
    # now it's time to parse the input file, which we have validated to check that it has
    # the columns we need. Each row will be parsed into an InputRow object
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        base_fieldnames = reader.fieldnames.copy()
        row_count = 1
        for row in reader:
            row_to_process = Genotype(row, args.amr_tool, organism_dict)
            # we only want to find matched rules for a row if it's relevant for AMR, so check this value first
            if row_to_process.to_process:                
                # extract the relevant rules for this ID, based on its organism
                relevant_rules = extract_relevant_rules(rules, row_to_process.organism)
                # determine if there's a matching rule for this row (this sets row_to_process.matched_rules)
                row_to_process.find_matching_rules(relevant_rules, amrfp_nodes)

            # track matched / unmatched hits for reporting
            # create a result row for each matched rule, as we need to duplicate rows in output
            # if they have multiple matching rules
            if row_to_process.matched_rules:
                matched_hits[row_count] = row_to_process.matched_rules
                for rule in row_to_process.matched_rules:
                    result_row = ResultRow(row_to_process, rule)
                    result_row.annotate_row(args.annot_opts)
                    result_row_objs.append(result_row)
            # if there's no matching rule, or it's a row we don't process, still create a ResultRow, but it has no rule
            else:
                unmatched_hits.append(row.get('Hierarchy node'))
                result_row = ResultRow(row_to_process, None)
                result_row.annotate_row(args.annot_opts)
                result_row_objs.append(result_row)

            # keep Genotype objects in case we need them later
            genotype_rows.append(row_to_process)
            row_count += 1
    

    # annotate each ResultRow and collect outputs for genotype report and summary
    #for rr in result_row_objs:
    #    base_row = dict(rr.base_row)
    #    base_row['organism'] = rr.organism
        # pass single-rule list or None to annotate_row
    #    annotated_rows = annotate_row(base_row, rr.rule_as_list(), args.annot_opts)
        # annotate_row returns a list; for a single-rule we expect one annotated dict
    #    if annotated_rows:
    #        # attach annotated dict to ResultRow for later summary processing
    #        rr.annotated_row = annotated_rows[0]
    #        output_rows.extend(annotated_rows)
    #        # build a summary-compatible dict via SummaryRow (supply conversion maps from ResourceManager)
    #        card_conv = rm().get_amrfp_card_conversion()
    #        card_map = rm().get_card_drug_class_map()
    #        sr = SummaryRow(rr, card_conversion=card_conv, card_drug_map=card_map)
    #        summary_input_rows.append(sr.to_summary_dict())
    #    else:
    #        # no annotated rows returned (shouldn't normally happen) -> still append base record
    #        rr.annotated_row = base_row
    #        output_rows.append(base_row)
    #        card_conv = rm().get_amrfp_card_conversion()
    #        card_map = rm().get_card_drug_class_map()
    #        sr = SummaryRow(rr, card_conversion=card_conv, card_drug_map=card_map)
    #        summary_input_rows.append(sr.to_summary_dict())
    
    # now write out the genotype report, which annotates each row with the rule info
    write_genotype_report(args, result_row_objs, unmatched_hits, matched_hits, base_fieldnames)
    
    # okay now we need the unique sample IDs, for the summary output
    # pass summary_input_rows (a list of summary-compatible dicts built from ResultRow/annotated rows)
    summary_output = prepare_summary(summary_input_rows, rules, samples_to_parse, args.no_flag_core, args.no_rule_interpretation)
    #summary_output = None

    write_genome_report(args, summary_output)


def download_resources():
    """
    Download and cache the AMRFP and CARD database files required.
    """
    rm().setup_all_resources()
    print("Resource download complete.")