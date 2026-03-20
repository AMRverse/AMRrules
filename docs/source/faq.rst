****
FAQ
****

.. dropdown:: I want to use AMRrules with input from <insert tool here>. Can I?

   Currently, AMRrules only supports input from AMRFinderPlus. We are working to include support for other tools and databases, see our :ref:`planned features <planned_other_tools>` for more information.

.. dropdown:: Why can't I find <insert rule> being applied in my output files?

   While the AMRrules engine itself currently only accepts AMRFinderPlus input, the rules themselves have been curated by organism experts who are aware of resistance mechanims in these organisms that may not yet be in the AMRFinderPlus database. Therefore, there are rules for some organisms which currently the engine cannot apply. Please also see the answer to the next question.

.. dropdown:: Why is there a rule for <insert mechanism> when AMRFinderPlus doesn't report it?

   We envision the rules as a place where experts can collate all currently knowledge of resistance mechanisms in their organism, and strongly support curators including rules for mechanisms which are currently not yet able to be detected by AMRFinderPlus. We work closely with the AMRFinderPlus team to help update their databases to include these new mechanisms, so future versions of the software can accurately report them.

.. dropdown:: Why isn't there a rule for <insert mechanism> in my favourite bug?

   Rules are regularly being updated and curated by the `ESGEM-AMR organism subgroups <https://esgem-amr.amrrules.org/>`__. If you have suggestions for new rules, please `post an issue on our GitHub <https://github.com/AMRverse/AMRrules/issues>`__, detailing your suggestion. We will pass this on to the relevant organism subgroup for consideration in the next update. Alternatively, you can also email us at esgem-amr [at] gmail [dot] com.

.. dropdown:: Why isn't my favourite bug supported?

   A full list of all organism subgroups can be found on the `ESGEM-AMR <https://esgem-amr.amrrules.org/>`__ website. If your organism of interest does not have a group, and you are interested in creating one, you can `register your interest here <https://docs.google.com/forms/d/e/1FAIpQLSeH96VlioxLKarZOLMqD-f1fLnb9WYOHYz4tZ9NtQzpHrKyzw/viewform?usp=sf_link>`__.

.. dropdown:: I think a rule is wrong. What do I do?

   If you think there is an error in a rule, we encourage you to first search our `GitHub issues page <https://github.com/AMRverse/AMRrules/issues>`__ to see if this rule is being discussed by others. If not, then please `post an issue <https://github.com/AMRverse/AMRrules/issues/new>`__ describing your concerns about the rule. Alternatively, you can also email us at esgem-amr [at] gmail [dot] com.

.. dropdown:: I don't know the species of my genome, what do I do?

   AMRrules interprets genotypes in an organism-specific manner, and so an organism is required for the engine to know what rules to apply. There are many methods you can use for determining the species of your genome. You could try uploading the genome to `PathogenWatch <https://pathogen.watch>`__ or using their `Speciator <https://github.com/pathogenwatch-oss/speciator>`__ tool to detect the species, or using tools provided by `GTDB <https://gtdb.ecogenomic.org/>`__.

.. dropdown:: I have a lot of partial hits in my result file, what does this mean?

   AMR genotyping and interpretation will only work as expected on good quality genome assemblies. The presence of lots of partial hits in your output file may suggest quality issues with your genome assembly, resulting in fragmentation of AMR gene hits across contigs. We recommend performing QC on your genome, you can find organism specific QC thresholds on `Qualibact <https://happykhan.github.io/qualibact/>`__.
