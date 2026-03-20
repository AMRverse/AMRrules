.. _planned_features:

Planned Features
================

This page outlines features and expansions that are planned for future releases of AMRrules.


Rules
~~~~~

* Expand rules to include all markers reported by `AMRFinderPlus <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>`_ for each supported organism. 
  
  * Where needed, we will proposed new markers for inclusion in the AMRFinderPlus database

* Rules will be developed to support both `EUCAST <https://www.eucast.org/>`_ and `CLSI <https://clsi.org/>`_ definitions of phenotypes, enabling organism-specific interpretation that aligns with established clinical breakpoint standards from both major guideline bodies.
 

Specification 
~~~~~~~~~~~~~

* We plan to expand the AMRrules Specification to include agreed quantitative fields, enabling more precise representation of resistance data beyond the current qualitative categories.

Engine
~~~~~~

.. _planned_other_tools:
* Extend compatibility to support additional widely-used AMR databases and tools, including:

  * `CARD <https://card.mcmaster.ca/>`_ (Comprehensive Antibiotic Resistance Database)
  * `ResFinder <https://cge.food.dtu.dk/services/ResFinder/>`_

* Some variation types are currently not by the AMRrules interpretation engine, as there are insufficient rules available for thorough testing. Support for the following variant types is planned for future development:

  * Gene truncation detected
  * Gene copy number variant detected
  * Nucleotide variant detected in multi-copy gene
  * Low frequency variant detected
