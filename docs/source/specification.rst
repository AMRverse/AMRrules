**************************
AMRrules Specification
**************************

.. _specification:
Rule Specification
=======

This section details how interpretive rules should be encoded in the AMRrules format.

The full list of fields is below, with guidelines on how each field should be specified and interpreted. 

This information is also available in a `Google sheet that includes the AMRrules template <https://docs.google.com/spreadsheets/d/1t6Lr_p-WAOY0yAXWKzoKk4yb56D2JdSqwImg4RZBvFA/edit?usp=sharing>`__, with allowed values encoded in drop-down menus, to facilitate rule curation. 

Full list of fields
^^^^^^^^^^^^^^

:download:`Download <tables/specifications_and_guidance.tsv>`

.. csv-table::
   :file: tables/specifications_and_guidance.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/specifications_and_guidance.tsv>`

AMRrules-specific controlled vocabularies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Variation type
--------------

Specifies the nature of the type of variation to which the rule applies. Based on the 'variant type' column in the hAMRonization specification, with additional terms from the NCIT ontology.

.. csv-table::
   :file: tables/variation_type_spec.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/variation_type_spec.tsv>`

Evidence codes
--------------

Specified using the `Evidence and Conclusion Ontology (ECO) <https://www.evidenceontology.org/>`__, this field indicates the nature of the evidence supporting the rule. More than one can be listed, and the field should include all forms of evidence available to support the rule (multiple entries separated with ', '). 

Any `ECO <https://www.evidenceontology.org/>`__ codes can be used, but curators are encouraged to choose from the subset listed here, which covers the types of evidence typically available to support resistance mechanisms in bacteria. Note the literature source for each type of evidence noted here should be indicated in the ``PMID`` field.

.. csv-table::
   :file: tables/evidence_codes.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/evidence_codes.tsv>`

Evidence grade
--------------

This fields indicates the expert curators’ overall assessment of the level of support provided by all evidence considered. It is modelled on the `GRADE <https://book.gradepro.org/>`__ (Grading of Recommendations, Assessment, Development, and Evaluation) approach to assessing the certainty of evidence to guide decision making in healthcare.

AMRrules aims to provide rules to interpret all markers that have been detected in a given species, but in many cases the evidence can be quite limited. The ``evidence grade`` field gives users an overall guide to the strength of evidence, and the ``evidence limitations`` field highlights what kind of evidence is lacking. 

Note that if no experimental evidence is available, the rule should NOT be graded as 'high', even if there is strong evidence of statistical association between genotype and phenotype in natural populations. (Future updates to the rule specification will include additional fields to record quantitative details of genotype/phenotype associations.)

There are four possible 'grades' in AMRrules, these are listed below with guidance on what they mean in the context of AMRrules (modelled on the `GRADE <https://gdt.gradepro.org/app/handbook/handbook.html#h.xr5ac2p2khuq>`__ framework).

.. csv-table::
   :file: tables/evidence_grades.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/evidence_grades.tsv>`


Evidence limitations
--------------

This fields indicates highlights what kind of evidence is lacking to support interpretation of this marker in this organism. All rules with an ``evidence grade`` other than 'high' should have at least one limitation recorded.

.. csv-table::
   :file: tables/evidence_limitations.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/evidence_limitations.tsv>`

Breakpoint condition
--------------

EUCAST, CLSI and others sometimes assign differen breakpoints for different clinical conditions, infection sites, or drug delivery routes (e.g. intravenous vs oral). In such cases, this field is used to indicate which specific breakpoint the rule was defined against. This will often be blank, indicating that the rule is not specific to any particular type of infection or delivery route.

The list of allowed terms is taken from the EUCAST and CLSI 2025 Breakpoints, sourced from the digitized versions in the `AMR R package <https://amr-for-r.org/>`__ using this command:
```
clinical_breakpoints %>% filter(guideline=="CLSI 2025" | guideline=="EUCAST 2025") %>% group_by(site) %>% count()
```

.. csv-table::
   :file: tables/breakpoint_condition_values.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/breakpoint_condition_values.tsv>`

Organism code
-------------

Each rule is assigned a ruleID, which starts with a 3-letter code to indicate the organism subgroup who curated the rule. The list of available organism subgroup codes is below.

.. csv-table::
   :file: tables/organism_subgroup_codes.tsv
   :delim: tab
   :header-rows: 1

:download:`Download <tables/organism_subgroup_codes.tsv>`

.. _syntax:
Variant Specification
=======

The AMRrules specification needs to be able to encode interpretive rules for all types of genetic variants relevant to AMR in bacteria.

In 2024, the `ESGEM-AMR <https://github.com/AMRverse/ESGEM-AMR>`__ working group collated and reviewed examples of known variants across diverse bacteria, and identified the following types of AMR variants:

* Gene presence detected
* Amino acid substitution or insertion
* Nucleotide substitution or insertion
* Gene truncated (loss of function)
* Mutation in promoter region (substitution, deletion or insertion, including IS)
* Gene copy number changes
* Mutations in multi-copy genes (e.g. 23S rRNA)
* Low frequency variants (i.e. heterozygosity)

It was concluded that all such variants could be adequately addressed using a combination of three fields:

* ``gene``
* ``mutation`` (based on `HGVS syntax <https://hgvs-nomenclature.org/stable/recommendations/summary/>`__, with some modifications)
* ``variation type`` (based on `hAMRonization <https://github.com/pha4ge/hAMRonization/tree/master/schema>`__ field `Genetic Variation Type <https://github.com/pha4ge/hAMRonization/blob/master/hAMRonization/constants.py>`__, with some additions).

Specific examples of each AMR variant are shown below, with proposed mutation syntax and variation types for each (note that other fields required for rule definition, like organism, refseq accession, context, PMID are not included here for simplicity, as they are not essential to illustrate how to define a specific *kind* of variation):

.. |dash| replace:: ``-``

========== ============ ============================= ================================================= ========================== ========== 
  ID         gene         mutation                      variation type                                    drug                       category           
========== ============ ============================= ================================================= ========================== ==========
  KPN0001    blaSHV       |dash|                        Gene presence detected                            ampicillin                 wt R               
  KPN0002    gyrA         p.Ser83Tyr                    Protein variant detected                          ciprofloxacin              nwt I              
  KPN0003    parC         p.Ser80Ile                    Protein variant detected                          ciprofloxacin              nwt I              
  KPN0004    ompK36       c.25C>T                       Nucleotide variant detected                       meropenem                  nwt S              
  KPN0005    ompK36       p.114_115insGlyAsp            Protein variant detected                          meropenem                  nwt I              
  KPN0006    mgrB         p.(1_100)                     Gene truncation detected                          colistin                   nwt R              
  ECO0001    ampC         c.-11C>T                      Promoter variant detected                         ceftriaxone                nwt R              
  ECO0002    ampC         c.-14_-13insGT                Promoter variant detected                         ceftriaxone                nwt R              
  ACI0001    blaOXA-58    c.(-35_1)ins[ISAba125:inv]    Promoter variant detected                         ceftriaxone                nwt R              
  NGO0002    23S rDNA     c.[2045A>G][3]                Nucleotide variant detected in multi-copy gene    azithromycin               nwt R              
  ECO0003    blaTEM       c.[3]                         Gene copy number variant detected                 piperacillin+tazobactam    nwt R              
  MTC0001    gyrA         p.[Ala94Gly][0.13]            Low frequency variant detected                    ciprofloxacin              nwt R              
========== ============ ============================= ================================================= ========================== ==========


Syntax for mutations
^^^^^^^^^^^^^^^^^^^^

Syntax for 'mutation' column follows `HGVS <https://hgvs-nomenclature.org/stable/>`__, including:

* Gene and protein start sites are position 1 (there is no position 0)
* Ranges are specified using ``x_y``; for insertions the coordinates are specified as inclusive_exclusive, otherwise ranges are inclusive_inclusive
* Unknown ranges are specified with parentheses, ``(x_y)``. E.g. ``p.(1_100)insGlyAsp`` means an insertion of 2 amino acids (Gly and Asp) anywhere between codons 1 and 100 inclusive (as opposed to a replacement of amino acids 1 through 100 with GlyAsp, which would be expressed as ``p.1_100delinsGlyAsp``).

  a. Coordinates are specified relative to the reference sequence of a protein (p) or coding sequence (c)

* Coordinates upstream of coding sequence are specified relative to the start site, with a hyphen, e.g. ``c.-35`` indicates 35 bp upstream
* Mutations in protein and DNA are specified differently, e.g.

  a. ``p.Ser83Tyr``: change to protein sequence from Ser to Tyr at codon 83
  b. ``c.25C>T``: change to nucleotide coding region from C to T at nucleotide position 25

* Stop codons are specified (in both DNA and protein variants) as ``Ter``
* Following `IUPAC <https://hgvs-nomenclature.org/stable/background/standards/#aacode>`__, ``X`` signifies any amino acid, ``N`` signifies any DNA base
* ``^`` (caret) is used as "or", e.g. ``p.(Gly719Ala^Ser)``
* The letters ``inv`` indicate the inverse (i.e. reverse complement) of a sequence
* Repeat sequences are specified as ``sequence[N]`` where ``N`` is the number of copies of the repeat

Syntax for ‘mutation’ specific to AMRrules
------------------------------------------

* AMRrules requires amino acids be specified as three-letter codes (whereas HGVS allows single-letter or three-letter codes)

  * Accordingly, the STOP codon should be specified as 'Ter' rather than '*'

* In HGVS you must specify the reference sequence explicitly using a sequence accession, followed by `:` and then the mutation, e.g. ``NF000285.3:p.Gly238Ser``. In AMRrules the gene is specified in separate column/s (‘gene’, ‘refseq accession’, ‘ARO accession’) and should not be repeated in the mutation column. So the above rule should be coded as:

  * gene = ``blaSHV``
  * node = ``blaSHV``
  * refseq accession = ``NF000285.3``
  * ARO accession = ``ARO:3000015``
  * mutation = ``p.Gly238Ser``

* In AMRrules, insertion sequences (IS) should be labelled with their IS name as per `ISfinder <https://isfinder.biotoul.fr/list_names_attributed.php>`__, as many do not have their own sequence accessions in refseq. E.g. insertion of ISAba125 should be specified as ``ins[ISAba125]``, and insertion in reverse orientation to the gene to which the rule applies should be specified as ``ins[ISAba125:inv]``.
* In AMRrules, rules intended to apply when a gene is present in a minimum of N copies can be specified using the ``[N]`` syntax to indicate the minimum repeat/copy number of the whole coding sequence, as ``c.[N]``. 

  a. Note this syntax does not convey any information about the location of the copies, i.e. ``c.[2]`` simply indicates that there are at least 2 copies of the gene detected in the genome, whether they are tandem repeats or in different replicons such as one in the chromosome and one in a plasmid.

* In HGVS, the presence of multiple alleles (i.e. heterozygous) is specified as a colon-separated list of allelic variants e.g. ``[allele1];[allele2]``. 
* In AMRrules, rules that apply to variation in a multi-copy gene can be specified in this way, with each allele explicitly stated.

  a. Alternatively if the rule applies when a minimum of N copies of the gene carry the mutation (e.g. mutation in ≥3 copies of 23S rRNA resulting in resistance to azithromycin), this can be abbreviated using the ``[N]`` syntax to indicate the minimum repeat/copy number, as ``c.[allele][N]`` or ``p.[allele][N]``, e.g. ``c.[2045A>G][3]``.
* In AMRrules, rules that apply to 'low frequency variants', i.e. when a minimum fraction of reads, P, support presence of the allelic variant in a sequenced population, the minimum fraction can be specified by extension of the syntax for copy number, as ``[X]``. E.g. ``p.[Ala94Gly][0.13]`` (`example <https://www.atsjournals.org/doi/full/10.1164/rccm.201703-0556OC>`__ from the *Mycobacterium tuberculosis gyrA* gene).

  a. To put another way, in AMRrules the repeat syntax ``[X]`` is interpreted as a minimum copy number if ``X`` is an integer, and as a minimum read fraction if ``X`` is a double/float between 0 and 1. 


Explanation of ‘mutation’ syntax relevant to known AMR variants
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``p.Ser83Tyr``: change to protein sequence from Ser to Tyr at codon 83
* ``c.25C>T``: change to nucleotide coding region from C to T at nucleotide position 25
* ``p.114_115insGlyAsp``: change to protein sequence, with an insertion of amino acids Gly and Asp between codons 114 and 115
* ``p.(1_100)``: truncation (of any kind) anywhere in the first 100 amino acids of the protein sequence
* ``c.-11C>T``: change to nucleotide sequence from C to T, 11 bases upstream of the start site for the gene.
* ``c.-14_-13insGT``: insertion of nucleotides GT between positions -14 and -13, upstream of the start site of the gene
* ``c.(-35_1)ins[ISAba125:inv]``: insertion of ISAba125, in reverse orientation (:inv), anywhere between 35 bases upstream of the start site, and the start of the gene coding sequence
* ``c.[2045A>G][3]``: substitution of A to G at position 2045 of the gene. This mutation must occur in minimum 3 copies
* ``c.[3]``: gene needs to be present with a minimum of 2 copies
* ``p.[Ala94Gly][0.13]``: protein variant is present in >13% of reads 


Combinatorial rules
^^^^^^^^^^^^^^^^^^^

Combinatorial rules are defined using logical expressions in the ‘gene’ column, where the objects of the expression are rule identifiers (``ruleID``) that can be used as shorthand labels for the variants defined by ``gene:mutation (variant type)`` specified in the corresponding rules. The ``variation type`` should be specified as ‘Combination’.

- Each rule must have a unique ``ruleID``, assigned by the curating subgroup and prefixed with a 3-letter code that identifies the subgroup. 
- E.g. in the table below, ``KPN0008`` can be used in a logical expression in the ‘gene’ column to demarcate ``gyrA:p.Ser83Tyr``, and ``KPN0013`` can be used to demarcate ``qnr (Gene presence detected)``.
- So, the combination of these two variants can be specified as ``KPN0008 & KPN0013``, which expands to ``gyrA:p.Ser83Tyr & qnr (Gene presence detected)``.

Rules **must** be specified explicitly if the effect of the combination is NOT the same as the ‘most resistant’ (in terms of exceeding breakpoints, R > I > S; or deviation from wildtype, nonwildtype > wildtype) predicted category of the component rules. E.g. in the table below:

- The individual rules ``KPN0008`` and ``KPN0009`` solo each have expected category ‘nonwildtype I’, but in combination we expect ‘nonwildtype R’, so we need to specify the rule for the combination ``KPN0008 & KPN0009``.
- The expected category for genomes meeting rule ``KPN0002`` (i.e. carrying core gene oqxA, => wildtype S) in addition to rule ``KPN0008`` (i.e. with an acquired gyrA mutation, => nonwildtype I) is nonwildtype I. This is the same, not greater, than one of the component rules (``KPN0008``) so we do not need to specify the combination explicitly.

Note this means the combination must be specified explicitly if the combined effect is LESS resistant than the ‘most resistant’ component, e.g. in `this example from TB <https://pubmed.ncbi.nlm.nih.gov/34460306/>`__, deletion in one gene renders the resistance mutation in another gene irrelevant so the combination must be specified.

========== ================================ ============= =========================== ================ ========== 
  ID         gene                             mutation      variation type              drug             category           
========== ================================ ============= =========================== ================ ========== 
  KPN0002    oqxA                             |dash|        Gene presence detected      ciprofloxacin    wt S               
  KPN0008    gyrA                             p.Ser83Tyr    Protein variant detected    ciprofloxacin    nwt I              
  KPN0009    parC                             p.Ser80Ile    Protein variant detected    ciprofloxacin    nwt I              
  KPN0013    qnr                              |dash|        Gene presence detected      ciprofloxacin    nwt I              
  KPN0051    KPN0008 & KPN0009                |dash|        Combination                 ciprofloxacin    nwt R              
  KPN0052    (KPN0008 | KPN0009) & KPN0013    |dash|        Combination                 ciprofloxacin    nwt R              
========== ================================ ============= =========================== ================ ==========  


