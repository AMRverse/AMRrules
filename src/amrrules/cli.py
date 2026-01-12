import argparse, os
from amrrules import rules_engine, __version__
from amrrules.utils import get_supported_organisms

def main():

    # Get list of valid organism names
    supported_organisms = get_supported_organisms()

    parser = argparse.ArgumentParser(description="Interpretation engine for AMRrules.")
    parser.add_argument('--input', type=str, help='Path to the input file.')
    parser.add_argument('--output-prefix', type=str, help='Prefix name for the output files.')
    parser.add_argument('--output-dir', '-d', type=str, default=os.getcwd(), help='Output directory. Default is current working directory.')
    parser.add_argument('--sample-id', type=str, help="If interpreting a single genome, can optionally provide sample ID here. If no sample_id is provided, and the first column of the input file doesn't define a sample_id, then the default value will be 'sample'.", default=None)

    org_args = parser.add_mutually_exclusive_group()
    org_args.add_argument('--organism', '-o', type=str, help=f"Organism to interpret. Use --list-organisms to see all supported organisms.")
    org_args.add_argument('--organism-file', '-of', type=str, help='Path to the organism file. This file should have two columns: genome name in col1 (matching the sample name in the first col of the input file), and col2 is the organism name, which should be one of the supported organisms. File should be in tab-delimited format, with no header')
    org_args.add_argument('--list-organisms', action='store_true', help='List all supported organisms and exit.')
    #TODO: implement card and resfinder options, currently only amrfp is supported
    parser.add_argument('--amr-tool', '-t', type=str, default='amrfp', help='AMR tool used to detect genotypes: options are amrfp, rgi, resfinder. Currently only amrfp is supported.')
    #parser.add_argument('--hamronized', '-H', action='store_true', help='Input file has been hamronized')
    # TODO: implement this option to allow for selection of different AMRFP databases
    #parser.add_argument('--amrfp_db_version', type=str, default='latest', help='Version of the AMRFP database used. Default is latest. NOTE STILL TO BE IMPLEMENTED')
    parser.add_argument('--no-rule-interpretation', '-nr', type=str, default = 'nwtR', choices=['nwtR', 'nwtS'], help='How to interpret hits that do not match a rule. Options are: nwtR (default) - all nonwildtype hits with no matching rule are interpreted as resistant; nwtS - all nonwildtype hits with no matching rule are interpreted as susceptible.')
    parser.add_argument('--annot-opts', '-a', type=str, default='minimal', choices=['minimal', 'full'], help='Annotation options: minimal (context, drug, phenotype, category, evidence grade), full (everything including breakpoints, standards, etc)')
    parser.add_argument('--flag-core', action='store_true', help='Turn on flagging core genes in the summary output')
    parser.add_argument('--full-disrupt', action='store_true', help='Show the full mutation detected by AMRFinderPlus for POINT_DISRUPT calls in the summary report, rather than just labelling them as gene:-')
    parser.add_argument('--print-non-amr', action='store_true', help='Include non-AMR rows (eg VIRULENCE, STRESS) from the input file in the interpreted output. By default, these rows are skipped.')
    parser.add_argument('--download-resources', action='store_true', help='Download AMRFinderPlus resource files and exit.')
    parser.add_argument('--version', action='version', version=f"amrrules {__version__}")

    args = parser.parse_args()

    if args.download_resources:
        rules_engine.download_resources()
        return
    
    if args.list_organisms:
        print("Supported organisms:")
        for org in supported_organisms:
            print(f"- {org}")
        return

    # Require input/output/organism for normal run
    if not args.input or not args.output_prefix or (not args.organism and not args.organism_file):
        parser.error('You must specify --input, --output_prefix, and --organism (or --organism_file) unless using --download-resources.')

    if args.amr_tool != 'amrfp':
        raise NotImplementedError("Currently only amrfp is supported. Please use amrfp as the AMR tool.")
    
    # can't have both organism_file and sample_id, as sample_id is only allowed for single entry files
    if args.organism_file and args.sample_id:
        parser.error("Please provide either --sample_id or --organism_file, not both. --sample_id should only be used if there is a single sample in the input file. Providing --organism_file presumes multiple samples to be processed, and is incompatible.")
    
    # check that the organism provided actually exists in the ruleset
    if args.organism:
        if args.organism not in supported_organisms:
            parser.error(f"Invalid organism name. Must be one of:\n{'\n'.join(supported_organisms)}")


    rules_engine.run(args)
