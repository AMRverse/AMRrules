.. toctree::
    :hidden:

   installation
   rules
   usage
   interpretation
   specification
   tests

**************************
AMRrules
**************************

About
=====

Organism-specific interpretation of antimicrobial susceptibility testing (AST) data is standard in clinical microbiology, with rules regularly reviewed by expert committees of `EUCAST <https://www.eucast.org/>`__ and `CLSI <https://clsi.org/>`__. We aim to provide an analagous resource to support organism-specific interpretation of antimicrobial resistance (AMR) genotypes derived from pathogen whole genome sequence (WGS) data.

AMRrules encode organism-specific rules for the interpretation of AMR genotype data. Users first take their genomes and run them through an AMR genotyping tool such as `AMRfinderplus <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>`__, which identifies AMR genes and mutations present in the genome. The resulting output can then interpreted using AMRrules, which applies the organism-specific rules to generate a final interpreted report of predicted antimicrobial susceptibility for that genome.

.. image:: images/amrfinder_pipeline.png
   :alt: Workflow for interpreting AMR genotypes using AMRrules
   :align: center

The AMRrules Python package includes the rules themselves (see `rules/` directory) as well as code to apply the rules to interpret AMR genotypes (currently limited to AMRFinderPlus output), generating informative genome reports that capture expert knowledge about how core and acquired genes and mutations contribute to antimicrobial susceptibility. 


Rule curation and development status
====================================

Rules are curated by organism experts belonging to `ESGEM-AMR <https://github.com/AMRverse/AMRrulesCuration/>`__, a working group of `ESGEM, the ESCMID Study Group on Epidemiological Markers <https://www.escmid.org/esgem/>`__. Initial rule curation has focused on defining rules for the interpretation of core genes and expected resistances, but acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.

Current rulesets and supported organisms can be found :ref:`here <rules>`.

The rule specification is available in `here <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`__ (v0.6, guidance on tab 2).

We are focusing early development on compatibility with NCBI resources (i.e. the `AMRfinderplus <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>`__ genotyping tool, and the associated NCBI databases including `AMR refgene <https://www.ncbi.nlm.nih.gov/pathogens/refgene/>`__, `AMR Reference Gene Hierarchy <https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy>`__, and the `Reference HMM Catalog <https://www.ncbi.nlm.nih.gov/pathogens/hmm/>`__). In future we plan for interoperability with `CARD <https://card.mcmaster.ca/>`__ and `ResFinder <http://genepi.food.dtu.dk/resfinder>`__ (and other tools based on these), using `hAMRonization <https://github.com/pha4ge/hAMRonization>`__.


Citation
========

TBD

Contributors
============
The AMRrules concept was initially workshopped by members of the `Holt lab <https://holtlab.net>`_ at `London School of Hygiene and Tropical Medicine <https://www.lshtm.ac.uk>`_ and further developed in collaboration with `Jane Hawkey <https://github.com/jhawkey>`_ at `Monash University <https://research.monash.edu/en/persons/jane-hawkey>`_. The AMRrules specification was developed by the ESGEM-AMR Data & Tools group, and the rules curated by the ESGEM-AMR Working Group (see `list of members <https://github.com/AMRverse/ESGEM-AMR/>`_), co-chaired by Natacha Couto (ESGEM Chair). 

Code was developed by Jane Hawkey and Kat Holt.
