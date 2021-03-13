Motivation
==========

The motivation for this project has its roots in Type 1 Diabetes (T1D) research. 
Delong Lab at CU Anshutz has discovered that Hybrid Insulin Peptides (HIPS) are a cause of T1D:
here_ is one of many papers regarding HIPS. Given the potential significance 
in being able to quickly identify hybrid peptides for this type of research, 
hybrid peptide identification tools are needed. As of now (December 2020), 
one other tools exists for finding hybrid peptides, NeoFusion_.

NeoFusion was built on top of MetaMorpheus. The authors of NeoFusion also 
wrote MetaMorpheus. Their goal was to leverage the power of MetaMorpheus to 
identify hybrid peptides. However, in their paper, they report not being able 
to identify trans-spliced peptides and only find cis-spliced peptides. In their
own analysis, roughly 20% of cis-spliced hybrid peptides were found using 
NeoFusion. Additionally, it takes 11 hours for their tool to run. Due to these 
limitations, there is both the space and the need for innovation.

.. _HIPS: https://pubmed.ncbi.nlm.nih.gov/30585061/
.. _NeoFusion: https://pubmed.ncbi.nlm.nih.gov/30346791/
