**************************
Usage
**************************

.. _usage:
Quickstart
==========

Interpreting AMRFinderPlus results for a single genome:

.. code-block:: bash

    # interpret AMRFinderPlus results for a single E. coli genome
    amrrules --input tests/data/input/test_ecoli_wildtype.tsv --output-prefix test_ecoli_wildtype --organism 's__Escherichia coli'

    # interpret AMRFinderPlus results for a single K. pneumoniae genome 
    amrrules --input tests/data/input/test_kpneumo_wildtype.tsv --output-prefix test_kpneumo_wildtype --organism 's__Klebsiella pneumoniae'


Interpret AMRFinderPlus results for multiple genomes of different organisms, by providing a sample/organism mapping file:

.. code-block:: bash

    amrrules  --input tests/data/input/test_multispp_amrfp.tsv --output-prefix test_multispp --organism-file tests/data/input/test_multispp_species.tsv


Usage
===============

Creating your AMRFinderPlus input
---------------------------------

AMRFinderPlus needs to be run with the ``--print_node`` option to include the necessary fields for AMRrules to match against. We recommend that you also provide the ``--name`` flag to ensure that AMRFinderPlus includes a column with the name of the sample you are running.

.. code-block:: bash

    amrfinder -n Kpn1.fasta --plus --print_node --name Kpn1 --organism Klebsiella_pneumoniae > Kpn1_AMRfp.tsv

Interpreting your AMRFinderPlus result
--------------------------------------

You can now interpret this output with AMRrules, using the *Klebsiella pneumoniae* rules by providing the organism name with the ``--organism`` flag. The organism name must be provided in the format ``s__Genus species``, and should be supplied in quotes. 

.. code-block:: bash

    amrrules --input Kpn1_AMRfp.tsv --output-prefix Kpn1_report --organism 's__Klebsiella pneumoniae'

.. note::

    Ensure you select the correct organism option here! We do not recommend applying rules for another organism, even if closely related to your target organism, as the rules have been curated in an organism-specific manner. 

In this example the user has also specified ``--plus`` when running AMRFinderPlus, which will return additional virulence, stress, or metal resistance genes detected by AMRFinderPlus. By default, these entries **will be omitted** from the final AMRrules interpreted output. If the user wishes to include these entries in the final output, add the ``--print-non-amr`` flag to the AMRrules commmand:

.. code-block:: bash

    amrrules --input Kpn1_AMRfp.tsv --output-prefix Kpn1_report --organism 's__Klebsiella pneumoniae' --print-non-amr

Naming samples in the output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We recommend users run AMRFinderPlus with the ``--name`` option to specify the sample name. This is particularly useful if you plan to interpret multiple samples at once using a single ``amrrules`` command. If you forget to set ``--name`` when running AMRFinderPlus, and have a single genome to interpret, you can provide the sample name to AMRrules using the ``--sample-id`` flag:

.. code-block:: bash

    amrrules --input Kpn1_AMRfp.tsv --output-prefix Kpn1_report --organism 's__Klebsiella pneumoniae' --sample-id Kpn1


Detailed options
=================

**Input & Output Options**

.. option:: --input <INPUT>

   Path to the tabular input file (must be AMRFinderPlus output for this version). Can be gzipped.

.. option:: --output-prefix <OUTPUT_PREFIX>

   Prefix name for generated output files.

.. option:: -d <OUTPUT_DIR>, --output-dir <OUTPUT_DIR>

   Output directory. **Default:** Current working directory.

.. option:: --sample-id <SAMPLE_ID>

   Sample identifier for single-genome interpretation. If omitted and the first column of the input file lacks a sample ID, defaults to ``sample``.

**Organism Selection**

.. option:: -o <ORGANISM>, --organism <ORGANISM>

   Organism name to interpret (e.g., `'s__Escherichia coli'`). Use :option:`--list-organisms` to see all supported targets.

.. option:: -of <ORGANISM_FILE>, --organism-file <ORGANISM_FILE>

   Path to a two-column, tab-delimited file without headers (`<sample_name>\t<organism_name>`).

.. option:: --list-organisms

   List all supported organisms and exit.

**Interpretation Settings**

.. option:: -t <AMR_TOOL>, --amr-tool <AMR_TOOL>

   AMR detection tool used to generate inputs. **Options:** ``amrfp`` (default).

.. option:: -nr {none,nwt,nwtR,nwtS}, --no-rule-interpretation {none,nwt,nwtR,nwtS}

   Handling for hits that do not match a rule. **Default:** ``none``. See :ref:`handling_ruleless_calls` for more details.

   * ``none``: No phenotype or clinical category assigned.
   * ``nwt``: Flagged as non-wildtype, but no clinical category set.
   * ``nwtR``: Interpreted as non-wildtype and marked Resistant.
   * ``nwtS``: Interpreted as non-wildtype and marked Susceptible.

.. option:: --trust-multi-copy

   Trust multi-copy gene calls without flagging them as potential draft assembly artifacts.

**Output Customization**

.. option:: -a {minimal,full}, --annot-opts {minimal,full}

   * ``minimal``: Context, drug, phenotype, category, and evidence grade.
   * ``full``: Includes breakpoints, standards, and extended metadata.

.. option:: --flag-core

   Flag core genes in the summary report output.

.. option:: --list-all-rules

   List all rules evaluated for a drug/class, regardless of whether they dictate the final result.

.. option:: --full-disrupt

   Report full mutation details for ``POINT_DISRUPT`` calls instead of truncating to ``gene:-``.

.. option:: --print-non-amr

   Include non-AMR targets (e.g., VIRULENCE, STRESS) in the output.

**Database & Resource Management**

.. option:: --download-resources [RESOURCE ...]

   Download or update reference databases. **Choices:** ``rules``, ``external``, or ``all`` (default).

.. option:: --db-version <DB_VERSION>

   Specific rules database release to download (e.g., ``v1.1.0``). Defaults to latest.

.. option:: --resource-info

   Display installed resource paths and database versions, then exit.


**Other**

.. option:: -h, --help

   Show help message and exit.

.. option:: --version

   Show program version number and exit.
