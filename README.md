# Interpretive standards for AMR genotypes

<img src="docs/source/images/AMRrules_logo.png" width="200" align="left">

Organism-specific interpretation of antimicrobial susceptibility testing (AST) data is standard in clinical microbiology, with rules regularly reviewed by expert committees of [EUCAST](https://www.eucast.org/) and [CLSI](https://clsi.org/). AMRrules aims to provide an analagous resource to support organism-specific interpretation of antimicrobial resistance (AMR) genotypes derived from pathogen whole genome sequence (WGS) data, curated by experts in the [ESGEM-AMR Working Group](https://amrverse.github.io/ESGEM-AMR/).

Users first take their genomes and run them through an AMR genotyping tool such as [AMRFinderPlus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/), which identifies AMR genes and mutations present in the genome. The resulting output can then interpreted using the AMRrules Python package, which applies the organism-specific rules to annotate each genotype maker and generate a final interpreted report of predicted antimicrobial susceptibility for that genome.

This AMRrules repository includes the rules themselves (see `rules/` directory) as well as the Python code to apply the rules to interpret AMR genotypes (currently limited to AMRFinderPlus output), generating informative genome reports that capture expert knowledge about how core and acquired genes and mutations contribute to antimicrobial susceptibility.

**[Full docs can be found here](https://amrrules.readthedocs.io/en/genome_summary_report_dev/index.html)**

[DOI: 10.5281/zenodo.12724317](https://doi.org/10.5281/zenodo.12724317)
