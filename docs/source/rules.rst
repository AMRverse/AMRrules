**************************
Rules
**************************

.. image:: images/organism_specific_rules.png
   :alt: Abbreviated example of organism-specific rules table
   :align: center
   :class: with-shadow

The full rule specification can be found at: `AMRrules spec v0.6 <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`__. Note this includes several additional fields beyond those pictured above, including NCBI and CARD ARO accessions to uniquely identify genes; details of the breakpoints and standards used; evidence codes, grades and limitations; and a rule annotation note.

Rule curation is a work in progress, under active development by the `ESGEM-AMR <https://github.com/AMRverse/AMRrulesCuration/>`__ Working Group.

Available rules
===============

Currently available rule sets are in the `rules/` directory of this repository, named by organism. In this beta release they focus mainly on core genes and expected resistances, however acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.

* *Acinetobacter baumannii* :doc:`../../rules/Acinetobacter_baumannii.txt`
* *Enterobacter* :doc:`../../rules/Enterobacter.txt`
* *Enterococcus faecalis* :doc:`../../rules/Enterococcus_faecalis.txt`
* *Enterococcus faecium* :doc:`../../rules/Enterococcus_faecium.txt`
* *Escherichia coli* :doc:`../../rules/Escherichia_coli.txt`
* *Klebsiella oxytoca* :doc:`../../rules/Klebsiella_oxytoca.txt`
* *Klebsiella pneumoniae* :doc:`../../rules/Klebsiella_pneumoniae.txt`
* *Neisseria gonorrhoeae* :doc:`../../rules/Neisseria_gonorrhoeae.txt` (includes acquired resistances, based on analysis of geno-pheno data)
* *Pseudomonas aeruginosa* :doc:`../../rules/Pseudomonas_aeruginosa.txt`
* *Salmonella* :doc:`../../rules/Salmonella.txt`
* *Staphylococcus aureus* :doc:`../../rules/Staphylococcus_aureus.txt`
* *Yersinia* :doc:`../../rules/Yersinia.txt`