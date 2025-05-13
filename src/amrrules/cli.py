import argparse
from amrrules import rules_engine

def main():
    parser = argparse.ArgumentParser(description="Interpreter for AMRrules.")
    parser.add_argument('--input', required=True)
    parser.add_argument('--output_prefix', required=True)
    parser.add_argument('--output_dir', default='.')
    parser.add_argument('--rules', '-r', required=True)
    parser.add_argument('--amr_tool', '-t', default='amrfp')
    parser.add_argument('--hamronized', '-H', action='store_true')
    parser.add_argument('--amrfp_db_version', default='latest')
    parser.add_argument('--annot_opts', '-a', default='minimal')

    args = parser.parse_args()
    rules_engine.run(args)
