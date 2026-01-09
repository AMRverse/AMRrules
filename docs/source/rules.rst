**************************
Rules and supported organisms
**************************

Available rules
===============

Currently available rule sets are in the ``rules/`` directory of this repository, named by organism. In this beta release they focus mainly on core genes and expected resistances, however acquired genes and mutations are included for some organisms already and will be added to others as the necessary data to define them accurately is accumulated and curated by the ESGEM-AMR working group.

Note that multiple organisms can exist within a rules file, particularly if the rules file belongs to a genus, or the named species belongs to a species complex (eg *Klebsiella oxytoca*). A full list of supported organisms can be found below.

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
* s__Enterobacter hormaechei
* s__Enterobacter roggenkampii
* s__Enterococcus faecalis
* s__Enterococcus faecium
* s__Escherichia coli
* s__Klebsiella grimontii
* s__Klebsiella huaxiensis
* s__Klebsiella michiganensis
* s__Klebsiella oxytoca
* s__Klebsiella pasteurii
* s__Klebsiella pneumoniae
* s__Klebsiella spallanzanii
* s__Neisseria gonorrhoeae
* s__Pseudomonas aeruginosa
* s__Salmonella enterica
* s__Staphylococcus aureus
* s__Yersinia enterocolitica
* s__Yersinia pseudotuberculosis

Rule specification
==================

The full rule specification can be found at: `AMRrules spec v0.6 <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`__. An abbreviated rules file is shown below.

.. image:: images/organism_specific_rules.png
   :alt: Abbreviated example of organism-specific rules table
   :align: center
   :class: with-shadow

Note the full specification includes several additional fields beyond those pictured above, including NCBI and CARD ARO accessions to uniquely identify genes; details of the breakpoints and standards used; evidence codes, grades and limitations; and a rule annotation note.

Rule curation is a work in progress, under active development by the `ESGEM-AMR <https://github.com/AMRverse/AMRrulesCuration/>`__ Working Group.