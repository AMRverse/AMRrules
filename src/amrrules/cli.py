import argparse, os
from amrrules import rules_engine, __version__
from amrrules.utils import get_supported_organisms, create_rules_file_key

def main():

    # Get list of valid organism names
    supported_organisms = get_supported_organisms()
    # create the rules key file
    create_rules_file_key()

    parser = argparse.ArgumentParser(description="Interpreter for AMRrules.")
    parser.add_argument('--input', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output_prefix', type=str, required=True, help='Prefix name for the output files.')
    parser.add_argument('--output_dir', '-d', type=str, default=os.getcwd(), help='Output directory. Default is current working directory.')
    # TODO: Implement file that lists file name (col1) and orgnaism (col2) to be used for the rules
    parser.add_argument('--organism', '-o', type=str, required=False, help=f"Organism to interpret. Must be one of: {', '.join(supported_organisms)}")
    # TODO: Implement this option to allow for multiple organisms
    parser.add_argument('--organism_file', '-of', type=str, required=False, help='Path to the organism file. This file should have two columns: genome name in col1 (matching the sample name in the first col of the input file), and col2 is the organism name, which should be one of the supported organisms. File should be in tab-delimited format, with no header')
    #TODO: implement card and resfinder options, currently only amrfp is supported
    parser.add_argument('--amr_tool', '-t', type=str, default='amrfp', help='AMR tool used to detect genotypes: options are amrfp, card, resfinder. Currently only amrfp is supported.')
    parser.add_argument('--hamronized', '-H', action='store_true', help='Input file has been hamronized')
    # TODO: implement this option to allow for selection of different AMRFP databases
    #parser.add_argument('--amrfp_db_version', type=str, default='latest', help='Version of the AMRFP database used. Default is latest. NOTE STILL TO BE IMPLEMENTED')
    parser.add_argument('--annot_opts', '-a', type=str, default='minimal', choices=['minimal', 'full'], help='Annotation options: minimal (context, drug, phenotype, category, evidence grade), full (everything including breakpoints, standards, etc)')
    parser.add_argument('--no_flag_core', action='store_true', help='Turn off flagging core genes in the summary output')
    parser.add_argument('--version', action='version', version=f"amrrules {__version__}")

    args = parser.parse_args()
    rules_engine.run(args)
