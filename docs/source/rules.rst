**************************
Rules and supported organisms
**************************

Available rules
===============

Currently available rule sets are in the ``rules/`` directory of this repository, named by organism. In this beta release they focus mainly on core genes and expected resistances, however acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.

* `Acinetobacter baumannii <https://github.com/AMRverse/AMRrules/blob/main/rules/Acinetobacter_baumannii.txt>`__
* `Enterobacter <https://github.com/AMRverse/AMRrules/blob/main/rules/Enterobacter.txt>`__
* `Enterococcus faecalis <https://github.com/AMRverse/AMRrules/blob/main/rules/Enterococcus_faecalis.txt>`__
* `Enterococcus faecium <https://github.com/AMRverse/AMRrules/blob/main/rules/Enterococcus_faecium.txt>`__
* `Escherichia coli <https://github.com/AMRverse/AMRrules/blob/main/rules/Escherichia_coli.txt>`__
* `Klebsiella oxytoca <https://github.com/AMRverse/AMRrules/blob/main/rules/Klebsiella_oxytoca.txt>`__
* `Klebsiella pneumoniae <https://github.com/AMRverse/AMRrules/blob/main/rules/Klebsiella_pneumoniae.txt>`__
* `Neisseria gonorrhoeae <https://github.com/AMRverse/AMRrules/blob/main/rules/Neisseria_gonorrhoeae.txt>`__ (includes acquired resistances, based on analysis of geno-pheno data)
* `Pseudomonas aeruginosa <https://github.com/AMRverse/AMRrules/blob/main/rules/Pseudomonas_aeruginosa.txt>`__
* `Salmonella <https://github.com/AMRverse/AMRrules/blob/main/rules/Salmonella.txt>`__
* `Staphylococcus aureus <https://github.com/AMRverse/AMRrules/blob/main/rules/Staphylococcus_aureus.txt>`__
* `Yersinia <https://github.com/AMRverse/AMRrules/blob/main/rules/Yersinia.txt>`__

Supported organisms
===================

A full list of supported organisms can be found by running ``amrrules --list-organisms``. Currently supported organisms are:

* s__Acinetobacter baumannii


An abbreviated rules file:

.. image:: images/organism_specific_rules.png
   :alt: Abbreviated example of organism-specific rules table
   :align: center
   :class: with-shadow

The full rule specification can be found at: `AMRrules spec v0.6 <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`__. Note this includes several additional fields beyond those pictured above, including NCBI and CARD ARO accessions to uniquely identify genes; details of the breakpoints and standards used; evidence codes, grades and limitations; and a rule annotation note.

Rule curation is a work in progress, under active development by the `ESGEM-AMR <https://github.com/AMRverse/AMRrulesCuration/>`__ Working Group.