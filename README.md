[![DOI](https://zenodo.org/badge/788956719.svg)](https://zenodo.org/doi/10.5281/zenodo.12724317)

# Interpretive standards for AMR genotypes

<img src="AMRrules_logo.png" width="200" align="left">

Organism-specific interpretation of antimicrobial susceptibility testing (AST) data is standard in clinical microbiology, with rules regularly reviewed by expert committees of [EUCAST](https://www.eucast.org/) and [CLSI](https://clsi.org/). We aim to provide an analagous resource to support organism-specific interpretation of antimicrobial resistance (AMR) genotypes derived from pathogen whole genome sequence (WGS) data.

AMRrules encode organism-specific rules for the interpretation of AMR genotype data, and are curated by organism experts belonging to [ESGEM-AMR](https://github.com/interpretAMR/AMRrulesCuration/), a working group of [ESGEM, the ESCMID Study Group on Epidemiological Markers](https://www.escmid.org/esgem/). The rule specification is available in [this Google sheet](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing) (v0.5, guidance on tab 2).

This AMRrules Python package includes the rules themselves (see `rules/` directory) as well as code to apply the rules to interpret AMR genotypes (currently limited to [AMRfinderplus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) output), generating informative genome reports that capture expert knowledge about how core and acquired genes and mutations contribute to antimicrobial susceptibility. 

We are focusing early development on compatibility with NCBI resources (i.e. the [AMRfinderplus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) genotyping tool, and the associated NCBI databases including [AMR refgene](https://www.ncbi.nlm.nih.gov/pathogens/refgene/), [AMR Reference Gene Hierarchy](https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy), and the [Reference HMM Catalog](https://www.ncbi.nlm.nih.gov/pathogens/hmm/)). In future we plan for interoperability with [CARD](https://card.mcmaster.ca/) and [ResFinder](http://genepi.food.dtu.dk/resfinder) (and other tools based on these), using [hAMRonization](https://github.com/pha4ge/hAMRonization).

Initial rule curation has focused on defining rules for the interpretation of core genes and expected resistances, but acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.

## Organism-specific rules
![rules_table](organism_specific_rules.png?raw=true)

Full specification: [AMRrules spec v0.5](https://docs.google.com/spreadsheets/d/1F-J-_8Kyo3W0Oh6eDYyd0N8ahqVwiddM2112-Fg1gKc/edit?usp=sharing). Note this includes several additional fields beyond those pictured above, and detailed guidance tabs.

### Available rules

Rule curation is being actively worked on by the ESGEM-AMR working group.

Currently available rule sets are in the [rules/](rules/) directory of this repository, named by organism, and focus mainly on core genes and expected resistances:

* [Acinetobacter baumannii](rules/Acinetobacter_baumannii.txt)
* [Escherichia coli](Escherichia_coli.txt)
* [Klebsiella pneumoniae](Klebsiella_pneumoniae.txt)
* [Neisseria gonorrhoeae](Neisseria_gonorrhoeae.txt) (acquired resistances, based on analysis of geno-pheno data)
* [Pseudomonas aeruginosa](Pseudomonas_aeruginosa.txt)
* [Salmonella](Salmonella.txt)
* [Staphylococcus aureus](Staphylococcus_aureus.txt)
* [Yersinia](Yersinia.txt)


## Interpreting AMRfinderplus genotype results with AMRrules: TO BE UPDATED

<img src="amrfinder_pipeline.png" width="600">

### Package installation

```
TBD
```

Use AMRrules to interpret AMRfinderplus results for a single genome (E. coli test data)

```
amrrules  --input tests/data/input/test_ecoli_genome.tsv
          --output_prefix test_ecoli_genome
          --organism 'Escherichia coli'
```


Use AMRrules to interpret results for multiple genomes of different organisms (test data)

```
amrrules  --input tests/data/input/test_data_amrfp_multiSpp.tsv
          --output_prefix test_multispp
          --organism_file 'test_data_sppCalls.tsv'
```


Example command to run AMRfinder plus, using the --print_node option

```
amrfinder -n example_data_kleb/ERR257656.fasta --plus
                                               --print_node
                                               --name ERR257656
                                               --organism Klebsiella_pneumoniae
                                               > ERR257656_amrFinderPlusOutput.tsv
```


## Example interpreted genotype report: TO BE UPDATED

## Example generic genome report: TO BE UPDATED
![genome_report](genome_report.png?raw=true)

Example file (PDF): [genome_report.pdf](genome_report.pdf)

Example file (RTF): [genome_report.rtf](genome_report.rtf)

## Contributors
The AMRrules concept was initially workshopped by members of the [Holt lab](https://holtlab.net) at [London School of Hygiene and Tropical Medicine](https://www.lshtm.ac.uk) and further developed in collaboration with [Jane Hawkey](https://github.com/jhawkey) at [Monash University](https://research.monash.edu/en/persons/jane-hawkey). The AMRrules specification was developed by the ESGEM-AMR Data & Tools group, and the rules curated by the ESGEM-AMR Working Group (see [list of members](https://github.com/interpretAMR/AMRrulesCuration/)), chaired by Natacha Couto (ESGEM Chair). Code was developed by Jane Hawkey and Kat Holt.
