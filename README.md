[![DOI](https://zenodo.org/badge/788956719.svg)](https://zenodo.org/doi/10.5281/zenodo.12724317)

# Interpretive standards for AMR genotypes

<img src="docs/source/images/AMRrules_logo.png" width="200" align="left">

Organism-specific interpretation of antimicrobial susceptibility testing (AST) data is standard in clinical microbiology, with rules regularly reviewed by expert committees of [EUCAST](https://www.eucast.org/) and [CLSI](https://clsi.org/). We aim to provide an analagous resource to support organism-specific interpretation of antimicrobial resistance (AMR) genotypes derived from pathogen whole genome sequence (WGS) data.

The AMRrules Python package includes the rules themselves (see `rules/` directory) as well as code to apply the rules to interpret AMR genotypes (currently limited to AMRFinderPlus output), generating informative genome reports that capture expert knowledge about how core and acquired genes and mutations contribute to antimicrobial susceptibility.

**[Full docs can be found here](https://amrrules.readthedocs.io/en/genome_summary_report_dev/index.html)**

AMRrules encode organism-specific rules for the interpretation of AMR genotype data. Users first take their genomes and run them through an AMR genotyping tool such as [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/), which identifies AMR genes and mutations present in the genome. The resulting output can then interpreted using AMRrules, which applies the organism-specific rules to generate a final interpreted report of predicted antimicrobial susceptibility for that genome.


