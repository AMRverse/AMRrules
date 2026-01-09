**************************
Usage
**************************

Quickstart
==========

Interprting AMRFinderPlus results for a single genome

::
    amrrules --input tests/data/input/test_ecoli_genome.tsv --output_prefix test_ecoli_genome --organism 's__Escherichia coli'
    amrrules --input tests/data/input/test_kleb_SGH10.tsv --output_prefix test_kleb_SGH10 --organism 's__Klebsiella pneumoniae'


Interpret AMRFinderPlus results for multiple genomes of different organisms

::
    amrrules  --input tests/data/input/test_data_amrfp_multiSpp.tsv --output_prefix test_multispp --organism_file tests/data/input/test_data_sppCalls.tsv


Example command to run AMRfinder plus (note the --print_node option)
and interpret the output with AMRrules

```
amrfinder -n Kpn1.fasta --plus --print_node --name Kpn1 --organism Klebsiella_pneumoniae > Kpn1_AMRfp.tsv

amrrules --input Kpn1_AMRfp.tsv --output_prefix Kpn1 --organism 's__Klebsiella pneumoniae'
```