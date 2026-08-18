from amrrules.rules_io import extract_unknown_core_rules, parse_rules_file, extract_relevant_rules
from amrrules.summariser import create_summary_dict
from amrrules.utils import check_sample_ids, validate_amrfp_file, get_organisms, open_input
from amrrules.output import write_genotype_report, write_genome_report
from amrrules.copy_number import apply_copy_number_rules, apply_combination_rules
from amrrules.resources import ResourceManager as rm
from amrrules.genotype_parser import GenoResult, Genotype
import csv
from importlib import resources
from collections import defaultdict

def run(args):

    # extract all the rules relevant to the organisms we're processing
    if args.organism_file:
        print("\nLoading organism assignments...")
        organism_dict, skipped_samples = get_organisms(args.organism_file)
        samples_with_org = set(organism_dict.keys())
    else:
        organism_dict = {'': args.organism}
        samples_with_org = None
        skipped_samples = None
    
    if args.amr_tool == 'amrfp':
        try:
            # then we need to grab the refgene heirarchy direct from the ncbi website (get latest for now)
            #TODO: user specifies version of amrfp database they used, or we extract this from hamronized file
            print("\nLoading AMRFinderPlus reference data...")
            amrfp_nodes = rm().refseq_nodes()
            # check the input file has the Hierarchy node column, and if an organism file is included, that the first column is Name
            samples_to_parse = validate_amrfp_file(args.input, multi_entry=bool(args.organism_file))
            
            # get the AMRFP to CARD conversion mapping for later use - we only want to do this once
            # so do it here and pass this to where it's needed
            card_amrfp_conversion = rm().get_amrfp_card_conversion() # get the AMRFP to CARD conversion mapping
            card_drug_map = rm().get_card_drug_class_map() # get CARD drugs and their associated classes
        except FileNotFoundError as exc:
            missing = f"\nMissing file: {exc.filename}" if getattr(exc, "filename", None) else ""
            raise SystemExit(
                "Required resource files were not found.\n"
                "Please run: amrrules --download-resources\n"
                "Then rerun your original command."
                f"{missing}\n"
                f"Details: {exc}"
            ) from None
    
    # if the is a multi-entry file, we need to check that all our sampleIDs are in the organism file
    # will raise an error if any are missing
    # will raise a warning if there are samples in the org file but aren't in the input file
    if args.organism_file:
        check_sample_ids(samples_with_org, samples_to_parse, skipped_samples)

    # collate unique rule files required for the organisms we need to parse
    rule_files = set()
    # open the rules key file and get the organism name
    key_file_path = resources.files("amrrules.rules").joinpath("rule_key_file.tsv")
    with open(key_file_path, 'r') as key_file:
        for row in key_file:
            # split the row into the organism and rules file
            organism, rules_filename = row.strip().split('\t')
            # if it's an organism we're interested in, add its rules file (only if it's not already)
            if organism in set(organism_dict.values()):
                rule_files.add(rules_filename)

    # parse the rule files
    print("\nParsing rule files...")
    rules = parse_rules_file(rule_files)

    matched_hits = {}
    unmatched_hits = []
    genotype_rows = []
    # now it's time to parse the input file, which we have validated to check that it has
    # the columns we need. Each row will be parsed into an InputRow object
    print("\nMatching markers to rules...")
    with open_input(args.input) as f:
        reader = csv.DictReader(f, delimiter='\t')
        base_fieldnames = reader.fieldnames.copy()
        row_count = 1
        for row in reader:
            if args.sample_id:
                row_to_process = GenoResult(row, args.amr_tool, organism_dict, args.print_non_amr, args.full_disrupt, sample_name=args.sample_id)
            else:
                row_to_process = GenoResult(row, args.amr_tool, organism_dict, args.print_non_amr, args.full_disrupt)
            # if this row belongs to a sample we should skip, update the to_process and to_print attributes to False
            if skipped_samples and row_to_process.sample_name in skipped_samples:
                row_to_process.to_process = False
                row_to_process.print_row = False
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
        if g.print_row:
            genotype_output_rows.extend(g.annotated_row)

    # we now want to create one object per rule/AMRFP subclass, so that we can summarise by drug or drug class.
    genotype_objects = []
    for g in genotype_rows:
        if g.to_process:
            if g.matched_rules:
                # go through each rule and create a duplicated Genotype object
                # however, only switch on duplicated for the second and subsequent rules, not the first one
                for rule in g.matched_rules:
                    if rule == g.matched_rules[0]:
                        geno_obj = Genotype.from_result_row(g, card_map=card_drug_map, rule=rule, duplicated = False)
                    else:
                        geno_obj = Genotype.from_result_row(g, card_map=card_drug_map, rule=rule, duplicated = True)
                    genotype_objects.append(geno_obj)
            else:
                # extract the subclasses and split as needed, assign drugs and classes that way
                g_subclasses = g.amrfp_subclass.split('/')
                for subclass in g_subclasses:
                    geno_obj = Genotype.from_result_row(g, card_amrfp=card_amrfp_conversion, amrfp_subclass=subclass, no_rule_interp=args.no_rule_interpretation)
                    genotype_objects.append(geno_obj)

    # now we want to group all of these objects by sample ID (if we have multiple samples)
    # because we need to apply copy number and combo rules by genome
    # and summarise by genome
    grouped_by_sample = defaultdict(list)
    for geno_obj in genotype_objects:
        grouped_by_sample[geno_obj.sample_name].append(geno_obj)

    # look for gene copy-number scenarios by sample, BEFORE drug-level
    # grouping. Runs against the full set of rules (not pre-filtered by
    # drug), since a copy-number rule's drug may have no other matched
    # rule at all for this sample
    for sample_name in grouped_by_sample:
        grouped_by_sample[sample_name] = apply_copy_number_rules(grouped_by_sample[sample_name], rules, card_drug_map)
        grouped_by_sample[sample_name] = apply_combination_rules(grouped_by_sample[sample_name], rules, card_drug_map)
        # sort AFTER both passes, so newly-appended rows are included
        #grouped_by_sample[sample_name].sort(key=lambda g: g.gene_symbol or '')

    # reorder the dict itself by sample_name, now that every sample's list is final
    #grouped_by_sample = dict(sorted(grouped_by_sample.items()))

    # Add new rows to the interpreted output, now that both copy number and combo
    # rules have been applied
    for sample_name, geno_objs in grouped_by_sample.items():
        for g in geno_objs:
            if getattr(g, 'copy_number_row', False) or getattr(g, 'combo_rule_row', False):
                genotype_output_rows.append(g.build_rule_only_row(args.annot_opts))

    genotype_output_rows.sort(key=lambda row: (row.get("Name", ""), row.get("gene", "")))

    # now write out the interpreted genotype report, which annotates each row with the rule info
    genotype_output_file = write_genotype_report(args, genotype_output_rows, base_fieldnames)

    # get the unknown rules, if there are any for this organism
    unknown_rules = extract_unknown_core_rules(rules, card_drug_map)

    # create the summary report
    summary_entry_dict = create_summary_dict(grouped_by_sample, unknown_rules, args.flag_core, args.no_rule_interpretation)
    # write out the summary report
    summary_output_file = write_genome_report(summary_entry_dict, args.output_dir, args.output_prefix)

    # print summary stats block
    num_skipped = len(skipped_samples) if skipped_samples is not None else 0
    ruler = "\u2500" * 52
    print()
    print(ruler)
    print(f"  \033[1;38;2;255;140;0mRun summary\033[0m")
    print(f"  Samples processed : {len(grouped_by_sample)}")
    print(f"  Samples skipped   : {num_skipped}")
    print(f"  Markers matched   : {len(matched_hits)}")
    print(f"  Markers unmatched : {len(unmatched_hits)}")
    print()
    print(f"  \033[1;32mOutput files\033[0m")
    print(f"  Interpreted genotype report   : {genotype_output_file}")
    print(f"  Genome summary report         : {summary_output_file}")
    print(ruler)
    print("\nAMRrules complete.")


def download_resources():
    """
    Download and cache the AMRFP and CARD database files required.
    """
    rm().setup_all_resources()
    print("Resource download complete.")