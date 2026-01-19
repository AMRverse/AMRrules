**************************
Usage
**************************

Quickstart
==========

Interpreting AMRFinderPlus results for a single genome:

::

    amrrules --input tests/data/input/test_ecoli_genome.tsv --output-prefix test_ecoli_genome --organism 's__Escherichia coli'
    amrrules --input tests/data/input/test_kleb_SGH10.tsv --output-prefix test_kleb_SGH10 --organism 's__Klebsiella pneumoniae'


Interpret AMRFinderPlus results for multiple genomes of different organisms, by providing a sample/organism mapping file:

::

    amrrules  --input tests/data/input/test_data_amrfp_multiSpp.tsv --output-prefix test_multispp --organism-file tests/data/input/test_data_sppCalls.tsv


Usage
===============

Creating your AMRFinderPlus input
---------------------------------

AMRFinderPlus needs to be run with the ``--print_node`` option to include the necessary fields for AMRrules to match against.

Eg:

::

    amrfinder -n Kpn1.fasta --plus --print_node --name Kpn1 --organism Klebsiella_pneumoniae > Kpn1_AMRfp.tsv

You can then interpret this output with AMRrules, using *Klebsiella pneumoniae* rules. Note that the organism name must be provided in the format ``s__Genus species``, and should be supplied in quotes.

::

    amrrules --input Kpn1_AMRfp.tsv --output-prefix Kpn1_report --organism 's__Klebsiella pneumoniae'

In this example the user has also specified ``--plus`` when running AMRFinderPlus, which will return additional virulence, stress, or metal resistance genes detected by AMRFinderPlus. By default, these entries **will be omitted** from the final AMRrules interpreted output. If the user wishes to include these entries in the final output, add the ``--print-non-amr`` flag to the AMRrules commmand:

::

    amrrules --input Kpn1_AMRfp.tsv --output-prefix Kpn1_report --organism 's__Klebsiella pneumoniae' --print-non-amr

Naming samples in the output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users may also wish to run AMRFinderPlus with the ``--name`` option to specify the sample name. This is particularly useful if you plan to interpret multiple samples at once using a single ``amrrules`` command. If you forget to set ``--name`` when running AMRFinderPlus, and have a single genome to interpret, you can provide the sample name to AMRrules using the ``--sample-id`` option:

::

    amrrules --input Kpn1_AMRfp.tsv --output-prefix Kpn1_report --organism 's__Klebsiella pneumoniae' --sample-id Kpn1


Detailed options
=================

::

    -h, --help            show this help message and exit
  --input INPUT         Path to the input file.
  --output-prefix OUTPUT_PREFIX
                        Prefix for the output files.
  --output-dir, -d OUTPUT_DIR
                        Output directory. Default is current working directory.
  --sample-id SAMPLE_ID
                        If interpreting a single genome, can optionally provide sample ID here. If no sample_id is provided, and the first column of the input file doesn't define a
                        sample_id, then the default value will be 'sample'.
  --organism, -o ORGANISM
                        Organism to interpret. Must be one of: s__Enterobacter hormaechei, s__Enterobacter roggenkampii, s__Enterococcus faecalis, s__Enterococcus faecium,
                        s__Klebsiella grimontii, s__Klebsiella huaxiensis, s__Klebsiella michiganensis, s__Klebsiella oxytoca, s__Klebsiella pasteurii, s__Klebsiella spallanzanii,
                        s__Salmonella enterica
  --organism-file, -of ORGANISM_FILE
                        Path to the organism file. This file should have two columns: genome name in col1 (matching the sample name in the first col of the input file), and col2 is
                        the organism name, which should be one of the supported organisms. File should be in tab-delimited format, with no header
  --list-organisms     List all supported organisms and exit.
  --amr-tool, -t AMR_TOOL
                        AMR tool used to detect genotypes: options are amrfp, rgi, resfinder. Currently only amrfp is supported.
  --no-rule-interpretation, -nr {nwtR,nwtS}
                        How to interpret hits that do not match a rule. Options are: nwtR (default) - all nonwildtype hits with no matching rule are interpreted as resistant; nwtS -
                        all nonwildtype hits with no matching rule are interpreted as susceptible.
  --annot-opts, -a {minimal,full}
                        Annotation options: minimal (context, drug, phenotype, category, evidence grade), full (everything including breakpoints, standards, etc)
  --flag-core           Turn on flagging core genes in the summary output
  --full-disrupt        Show full POINT_DISRUPT mutation details in the summary output, rather than just indicating gene inactivation with '-'.
  --print-non-amr       Include non-AMR rows (eg VIRULENCE, STRESS) from the input file in the interpreted output. By default, these rows are skipped.
  --download-resources  Download AMRFinderPlus resource files and exit.
  --version             show AMRrules version number and exit
