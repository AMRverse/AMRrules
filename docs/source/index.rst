.. toctree::
   :hidden:

=================
AMRrules
=================

AMRrules encode organism-specific rules for the interpretation of AMR genotype data, and are curated by organism experts belonging to `ESGEM-AMR <https://github.com/AMRverse/AMRrulesCuration/>`_, a working group of `ESGEM, the ESCMID Study Group on Epidemiological Markers <https://www.escmid.org/esgem/>`_. The rule specification is available in `this Google sheet <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`_ (v0.6, guidance on tab 2).

Organism-specific interpretation of antimicrobial susceptibility testing (AST) data is standard in clinical microbiology, with rules regularly reviewed by expert committees of `EUCAST <https://www.eucast.org/>`_ and `CLSI <https://clsi.org/>`_. We aim to provide an analagous resource to support organism-specific interpretation of antimicrobial resistance (AMR) genotypes derived from pathogen whole genome sequence (WGS) data.

This AMRrules Python package includes the rules themselves (see `rules/` directory) as well as code to apply the rules to interpret AMR genotypes (currently limited to `AMRfinderplus <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>`_ output), generating informative genome reports that capture expert knowledge about how core and acquired genes and mutations contribute to antimicrobial susceptibility. 

We are focusing early development on compatibility with NCBI resources (i.e. the `AMRfinderplus <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>`_ genotyping tool, and the associated NCBI databases including `AMR refgene <https://www.ncbi.nlm.nih.gov/pathogens/refgene/>`_, `AMR Reference Gene Hierarchy <https://www.ncbi.nlm.nih.gov/pathogens/genehierarchy>`_, and the `Reference HMM Catalog <https://www.ncbi.nlm.nih.gov/pathogens/hmm/>`_). In future we plan for interoperability with `CARD <https://card.mcmaster.ca/>`_ and `ResFinder <http://genepi.food.dtu.dk/resfinder>`_ (and other tools based on these), using `hAMRonization <https://github.com/pha4ge/hAMRonization>`_.

Initial rule curation has focused on defining rules for the interpretation of core genes and expected resistances, but acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.