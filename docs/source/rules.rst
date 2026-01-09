**************************
Rules
**************************

.. image:: extras/organism_specific_rules.png
   :alt: Abbreviated example of organism-specific rules table
   :align: center
   :class: with-shadow

The full rule specification can be found at: `AMRrules spec v0.6 <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`__. Note this includes several additional fields beyond those pictured above, including NCBI and CARD ARO accessions to uniquely identify genes; details of the breakpoints and standards used; evidence codes, grades and limitations; and a rule annotation note.

Rule curation is a work in progress, under active development by the `ESGEM-AMR <https://github.com/AMRverse/AMRrulesCuration/>`__ Working Group.

Available rules
=============================

Currently available rule sets are in the `rules/` directory of this repository, named by organism. In this beta release they focus mainly on core genes and expected resistances, however acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.

* [Acinetobacter baumannii](rules/Acinetobacter_baumannii.txt)
* [Enterobacter](rules/Enterobacter.txt)
* [Enterococcus faecalis](rules/Enterococcus_faecalis.txt)
* [Enterococcus faecium](rules/Enterococcus_faecium.txt)
* [Escherichia coli](rules/Escherichia_coli.txt)
* [Klebsiella pneumoniae](rules/Klebsiella_pneumoniae.txt)
* [Neisseria gonorrhoeae](rules/Neisseria_gonorrhoeae.txt) (acquired resistances, based on analysis of geno-pheno data)
* [Pseudomonas aeruginosa](rules/Pseudomonas_aeruginosa.txt)
* [Salmonella](rules/Salmonella.txt)
* [Staphylococcus aureus](rules/Staphylococcus_aureus.txt)
* [Yersinia](rules/Yersinia.txt)