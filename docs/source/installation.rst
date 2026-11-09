**************************
Installation
**************************

Dependencies
=============

AMRrules requires: 

* Python >= 3.12
* obonet >= 1.1.1 (to parse the CARD ontology)


Download and install AMRrules
=============================

Currently, AMRrules is only available for installation via source. We recommend you set up a new conda environment.

.. code-block:: bash

    # create your conda environment
    conda create -n amrrules -c bioconda python=3.12 pip
    conda activate amrrules

    # clone the repository
    git clone https://github.com/AMRverse/AMRrules
    cd AMRrules

    # install AMRrules
    make dev

After installation, you must download the AMRrulesDB (contains all the necessary rule sets), as well as key CARD and AMRFinderPlus resource files. Run:

.. code-block:: bash

    amrrules --download-resources
    
This will download the relevant databases required for AMRrules to function. You only need to run this **once** after installation, or when updating rules or the external databases (e.g. a new AMRFinderPlus database has been released).

You can update either the rules or the external databases, or all, by running:

.. code-block:: bash

    # download everything (rules + external databases)
    amrrules --download-resources all
    
    # download latest version of the rules
    amrrules --download-resources rules
    
    # download a specific version of the rules (e.g. v1.0.0)
    amrrules --download-resources rules --db-version 1.0.0
    
    # download new versions of the external databases (CARD, AMRFinderPlus)
    amrrules --download-resources external

More info on usage options can be found on the :ref:`Usage <usage>` page.

.. _check-the-installation:

Check the installation
======================

Check which organisms have rule-sets available in the installation:

.. code-block:: bash

    amrrules --list-organisms

Test on one of the included example datasets:
    
.. code-block:: bash

    amrrules --input tests/data/input/test_ecoli_wildtype.tsv --output-prefix test_ecoli_wildtype --organism 's__Escherichia coli'

Compare these resulting output files to those in ``tests/data/example_output``.

**Interpreted file**

.. code-block:: bash
    
    tests/data/example_output/test_ecoli_wildtype_interpreted.tsv

.. container:: table-wrapper

    .. csv-table::
      :file: ../../tests/data/example_output/test_ecoli_wildtype_interpreted.tsv
      :delim: tab
      :header-rows: 1
      :class: full-width

**Genome summary file**

.. code-block:: bash

    tests/data/example_output/test_ecoli_wildtype_genome_summary.tsv

.. container:: table-wrapper

    .. csv-table::
        :file: ../../tests/data/example_output/test_ecoli_wildtype_genome_summary.tsv
        :delim: tab
        :header-rows: 1
        :class: full-width

Addtional test data
========================

Additional test data and example commands can be found on the :ref:`Tests <test>` page.


