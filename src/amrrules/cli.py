import argparse, os
from amrrules import rules_engine, __version__
from amrrules.utils import get_supported_organisms

def main():

    # Get list of valid organism names
    supported_organisms = get_supported_organisms()


    parser = argparse.ArgumentParser(description="Interpretation engine for AMRrules.")
    parser.add_argument('--input', type=str, help='Path to the input file.')
    parser.add_argument('--output_prefix', type=str, help='Prefix name for the output files.')
    parser.add_argument('--output_dir', '-d', type=str, default=os.getcwd(), help='Output directory. Default is current working directory.')
    parser.add_argument('--download-resources', action='store_true', help='Download AMRFinderPlus resource files and exit.')

    org_args = parser.add_mutually_exclusive_group()
    org_args.add_argument('--organism', '-o', type=str, help=f"Organism to interpret. Must be one of: {', '.join(supported_organisms)}")
    org_args.add_argument('--organism_file', '-of', type=str, help='Path to the organism file. This file should have two columns: genome name in col1 (matching the sample name in the first col of the input file), and col2 is the organism name, which should be one of the supported organisms. File should be in tab-delimited format, with no header')
    #TODO: implement card and resfinder options, currently only amrfp is supported
    parser.add_argument('--amr_tool', '-t', type=str, default='amrfp', help='AMR tool used to detect genotypes: options are amrfp, card, resfinder. Currently only amrfp is supported.')
    parser.add_argument('--hamronized', '-H', action='store_true', help='Input file has been hamronized')
    # TODO: implement this option to allow for selection of different AMRFP databases
    #parser.add_argument('--amrfp_db_version', type=str, default='latest', help='Version of the AMRFP database used. Default is latest. NOTE STILL TO BE IMPLEMENTED')
    parser.add_argument('--no-rule-interpretation', '-nr', type=str, default = 'nwtR', choices=['nwtR', 'nwtS'], help='How to interpret hits that do not match a rule. Options are: nwtR (default) - all nonwildtype hits with no matching rule are interpreted as resistant; nwtS - all nonwildtype hits with no matching rule are interpreted as susceptible.')
    parser.add_argument('--annot_opts', '-a', type=str, default='minimal', choices=['minimal', 'full'], help='Annotation options: minimal (context, drug, phenotype, category, evidence grade), full (everything including breakpoints, standards, etc)')
    parser.add_argument('--no_flag_core', action='store_true', help='Turn off flagging core genes in the summary output')
    parser.add_argument('--version', action='version', version=f"amrrules {__version__}")

    args = parser.parse_args()

    if args.download_resources:
        rules_engine.download_resources()
        return

    # Require input/output/organism for normal run
    if not args.input or not args.output_prefix or (not args.organism and not args.organism_file):
        parser.error('You must specify --input, --output_prefix, and --organism or --organism_file unless using --download-resources.')

    if args.amr_tool != 'amrfp':
        raise NotImplementedError("Currently only amrfp is supported. Please use amrfp as the AMR tool.")

    rules_engine.run(args)
