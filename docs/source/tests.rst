**************************
Test Data
**************************

Test installation
==========

To test that AMRrules is installed and working correctly, test data is included inside ``tests/data/input/``. This includes examples of AMRFinderPlus output files for *Escherichia coli* and *Klebsiella pneumoniae* genomes, as well as a multi-species example file.

To run the AMRrules test data, use the following commands:

::

    # test a single E. coli genome
    amrrules --input tests/data/input/test_ecoli_genome.tsv --output-prefix test_ecoli_genome --organism 's__Escherichia coli'
    
    # test a single wild-type K. pneumoniae genome, and flag core genes in the summary output
    amrrules --input tests/data/input/test_kleb_SGH10.tsv --output-prefix test_kleb_SGH10 --organism 's__Klebsiella pneumoniae' --flag-core

    # test a multi-genome, multi-species example
    amrrules  --input tests/data/input/test_data_amrfp_multiSpp.tsv --output-prefix test_multispp --organism-file tests/data/input/test_data_sppCalls.tsv


Outputs will be written to your current working directory, and can be compared to the output files in ``tests/data/example_output``.

Beta-testing data
=================


