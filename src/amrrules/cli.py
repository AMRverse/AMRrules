import argparse, os
from amrrules import rules_engine, __version__

def main():
    parser = argparse.ArgumentParser(description="Interpreter for AMRrules.")
    parser.add_argument('--input', type=str, required=True, help='Path to the input file.')
    parser.add_argument('--output_prefix', type=str, required=True, help='Prefix name for the output files.')
    parser.add_argument('--output_dir', '-d', type=str, default=os.getcwd(), help='Directory to write the output files to. Default is current directory.')
    # TODO: Implement file that lists file name (col1) and orgnaism (col2) to be used for the rules
    parser.add_argument('--organism', '-o', type=str, required=True, help='Organism name of rules to apply.')
    #parser.add_argument('--rules', '-r', type=str, required=True, help='Path to the rules file, in tab-delimited format.')
    #TODO: implement card and resfinder options, currently only amrfp is supported
    parser.add_argument('--amr_tool', '-t', type=str, default='amrfp', help='AMR tool used to detect genotypes: options are amrfp, card, resfinder. Currently only amrfp is supported.')
    parser.add_argument('--hamronized', '-H', action='store_true', help='Input file has been hamronized')
    parser.add_argument('--amrfp_db_version', type=str, default='latest', help='Version of the AMRFP database used. Default is latest. NOTE STILL TO BE IMPLEMENTED')
    parser.add_argument('--annot_opts', '-a', type=str, default='minimal', help='Annotation options: minimal (context, drug, phenotype, category, evidence grade), full (everything including breakpoints, standards, etc)')
    parser.add_argument('--version', action='version', version=f"amrrules {__version__}")

    args = parser.parse_args()
    rules_engine.run(args)
