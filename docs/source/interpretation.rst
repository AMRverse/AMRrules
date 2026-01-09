**************************
Interpreting genotypes
**************************

When given an AMRFinderPlus input file, and an organism, AMRrules will apply the most relevant organism-specific rule to each AMRFinderPlus hit.

When annotating a hit with a rule, AMRrules will add the following information from the matching rule to the AMRFinderPlus output:
* ruleID
* gene context
* drug
* drug class
* phenotype
* clinical category
* evidence grade
* version
* organism

To add the full rule annotation information, set ``--annot_opts full``, which will add the following additional fields:
* breakpoint
* breakpoint standard
* breakpoint condition
* evidence code
* evidence limitations
* PMID
* rule curation note


By default, if no organism-specific rule is found for a hit, AMRrules will apply a generic call based on the setting for ``--no-rule-interpretation``. By default, this is set to ``nwtR``: any unmatched hit will be given the phenotype ``nonwildtype`` and the clinical category of resistance will be set to ``R``. Users can change this to be ``nwtS``, if they do not wish to call unmatched hits as resistant.