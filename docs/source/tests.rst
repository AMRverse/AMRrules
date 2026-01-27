**************************
Test Data
**************************

.. _test:
Test installation
=================

To test that AMRrules is installed and working correctly, test data is included inside ``tests/data/input/``. This includes examples of AMRFinderPlus output files for *Escherichia coli* and *Klebsiella pneumoniae* genomes, as well as a multi-species example file.

To run the AMRrules test data, use the following commands:

::

    # test a single E. coli genome (wildtype strain with no acquired resistance phenotypes)
    amrrules --input tests/data/input/test_ecoli_wildtype.tsv --output-prefix test_ecoli_wildtype --organism 's__Escherichia coli'
    
    # test a single wild-type K. pneumoniae genome, and flag core genes in the summary output
    amrrules --input tests/data/input/test_kleb_wildtype.tsv --output-prefix test_kleb_wildtype --organism 's__Klebsiella pneumoniae' --flag-core

    # test a multi-genome, multi-species example
    amrrules  --input tests/data/input/test_multispp_amrfp.tsv --output-prefix test_multispp --organism-file tests/data/input/test_multispp_species.tsv


Outputs will be written to your current working directory, and can be compared to the output files in ``tests/data/example_output``.

Beta-testing data
=================

Additional test datasets are available for select organisms.

**E. coli**: 20 *E. coli* genomes ranging from wildtype susceptible to multidrug resistant.

Genomes sourced from Mills et al, 2022, which includes AST profiles available `here <https://doi.org/10.1186/s13073-022-01150-7>`__.

:: 

    amrrules --input tests/data/input/test_ecoli_20strains.tsv --output-prefix test_ecoli_20strains --organism 's__Escherichia coli'


