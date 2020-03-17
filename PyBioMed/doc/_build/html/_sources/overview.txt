..  -*- coding: utf-8 -*-

Overview
========

To develop a powerful model for prediction tasks by machine learning algorithms such as sckit-learn, one of the most important things to consider is how to effectively represent the molecules under investigation such as small molecules, proteins, DNA and even complex interactions, by a descriptor. PyBioMed is a feature-rich package used for the characterization of various complex biological molecules and interaction samples, such as chemicals, proteins, DNA, and their interactions. PyBioMed calculates nine types of features including chemical descriptors or molecular fingerprints, structural and physicochemical features of proteins and peptides from amino acid sequence, composition and physicochemical features of DNA from their primary sequences, chemical-chemical interaction features, chemical-protein interaction features, chemical-DNA interaction features, protein-protein interaction features, protein-DNA interaction features, and DNA-DNA interaction features. We hope that the package can be used for exploring questions concerning structures, functions and interactions of various molecular data in the context of chemoinformatics, bioinformatics, and systems biology.  The python package PyBioMed is designed by CBDD Group (Computational Biology & Drug Design Group), Xiangya School of Pharmaceutical Sciences, Central South University.

Who uses PyBioMed?
~~~~~~~~~~~~~~~~~

For those researchers from different biomedical fields, the PyBioMed package can be used to analyze and represent various complex molecular data under investigation. PyBioMed will be helpful when exploring questions concerning structures, functions and interactions of various molecular data in the context of chemoinformatics, bioinformatics, and systems biology.


Motivation
~~~~~~~~~~
PyBioMed is intended to provide

-  Tools for pretreating molecules, proteins sequence and DNA sequence

-  Calculating chemical descriptors or molecular fingerprints from
   molecules' structures

-  Calculating structural and physicochemical features of proteins and peptides
   from amino acid sequence

-  Calculating composition and physicochemical features of DNA
   from their primary sequences

-  Calculating interaction features including chemical-chemical interaction features,
   chemical-protein interaction features, chemical-DNA interaction features,
   protein-protein interaction features, protein-DNA interaction features
   and DNA-DNA interaction features.

-  Getting molecular structures, protein sequence and DNA sequence from Internet through
   the molecular ID, protein ID and DNA ID.

Feature overview
~~~~~~~~~~~~~~~~

The table below shows the descriptors and the number of the descriptor that PyBioMed can calculate in four modules including PyMolecule, PyProtein, PyDNA and PyInteraction. PyMolecule module can calculate 14 different types of molecular descriptors and 18 different types of molecular fingerprints. PyProtein module can calculate 14 types of protein descriptors. PyDNA module can calculate 14 types of DNA descriptors and the number in the table appears when parameters are 'all_property = True, lamada=2, w=0.05'. PyInteraction module can calculate three types of descriptors.

+------------------+-------------------------------------------------------+--------------+
|Types             |Features                                               |Description   |
+==================+=======================================================+==============+
|PyMolecule        | - Constitution (30)                                   |              |
|                  | - Connectivity descriptors (44)                       |              |
|                  | - Topology descriptors (35)                           |              |
|                  | - Basak descriptors (21)                              |              |
|                  | - Burden descriptors (64)                             |              |
|                  | - Kappa descriptors (7)                               |              |
|                  | - E-state descriptors (237)                           |              |
|                  | - Moran autocorrelation descriptors (32)              |              |
|                  | - Geary autocorrelation descriptors (32)              |              |
|                  | - Molecular property descriptors (6)                  |              |
|                  | - Moreau-Broto autocorrelation descriptors (32)       |              |
|                  | - Charge descriptors (25)                             |`PyMolecule`_ |
|                  | - MOE-type descriptors (60)                           |              |
|                  | - CATS2D descriptors (150)                            |              |
+                  +-------------------------------------------------------+              +
|                  | - Daylight-type fingerprints (2048)                   |              |
|                  | - MACCS fingerprints (166)                            |              |
|                  | - Atom pairs fingerprints (1024)                      |              |
|                  | - TopologicalTorsion fingerprints (1024)              |              |
|                  | - E-state fingerprints (79)                           |              |
|                  | - FP2 fingerprints (1024)                             |              |
|                  | - FP3 fingerprints (210)                              |              |
|                  | - FP4 fingerprints (307)                              |              |
|                  | - ECFP2 fingerprints (1024)                           |              |
|                  | - ECFP4 fingerprints (1024)                           |              |
|                  | - ECFP6 fingerprints  (1024)                          |              |
|                  | - Morgan fingerprints (1024)                          |              |
|                  | - Ghosecrippen fingerprints (110)                     |              |
|                  | - FCFP2 fingerprints (1024)                           |              |
|                  | - FCFP4 fingerprints (1024)                           |              |
|                  | - FCFP6 fingerprints (1024)                           |              |
|                  | - Pharm2D2point fingerprints (135)                    |              |
|                  | - Pharm2D3point fingerprints (2135)                   |              |
+------------------+-------------------------------------------------------+--------------+
|PyProtein         | - Amino acid composition (20)                         |              |
|                  | - Dipeptide composition (400)                         |              |
|                  | - Tripeptide composition (8000)                       |              |
|                  | - CTD composition (21)                                |              |
|                  | - CTD transition (21)                                 |              |
|                  | - CTD distribution (105)                              |              |
|                  | - M-B autocorrelation (240)                           |`PyProtein`_  |
|                  | - Moran autocorrelation (240)                         |              |
|                  | - Geary autocorrelation (240)                         |              |
|                  | - Conjoint triad features (343)                       |              |
|                  | - Quasi-sequence order descriptors (100)              |              |
|                  | - Sequence order coupling number (60)                 |              |
|                  | - Pseudo amino acid composition 1 (50)                |              |
|                  | - Pseudo amino acid composition 2 (50)                |              |
+------------------+-------------------------------------------------------+--------------+
|PyDNA             | - Basic kmer (16) (k=2)                               |              |
|                  | - Reverse compliment kmer (10) (k=2)                  |              |
|                  | - DAC (76) (all_property=True)                        |              |
|                  | - DCC (2812) (all_property=True)                      |              |
|                  | - DACC (2888) (all_property=True)                     |              |
|                  | - TAC (24) (all_property=True)                        |              |
|                  | - TCC (264) (all_property=True)                       |              |
|                  | - TACC (288) (all_property=True)                      |`PyDNA`_      |
|                  | - PseDNC (18) (all_property=True,lamada=2,w=0.05)     |              |
|                  | - PseKNC (66) (all_property=True,lamada=2,w=0.05)     |              |
|                  | - PC-PseDNC (18) (all_property=True,lamada=2,w=0.05)  |              |
|                  | - PC-PseTNC (66) (all_property=True,lamada=2,w=0.05)  |              |
|                  | - SC-PseDNC (92) (all_property=True,lamada=2,w=0.05)  |              |
|                  | - SC-PseTNC (88) (all_property=True,lamada=2,w=0.05)  |              |
+------------------+-------------------------------------------------------+--------------+
|PyInteraction     | - Feature type 1                                      |              |
|                  | - Feature type 2                                      |`PyInter`_    |
|                  | - Feature type 3                                      |              |
+------------------+-------------------------------------------------------+--------------+




The Python programming language
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Python is a powerful programming language that allows simple and flexible representations of biochemical molecules, and clear and concise expressions of bioinformatics algorithms. Python has a vibrant and growing ecosystem of packages that PyBioMed uses to provide more features such as RDkit and Pybel. In addition, Python is also an excellent “glue” language for putting together pieces of software from other languages which allows reuse of legacy code and engineering of high-performance algorithms. Equally important, Python is free, well-supported, and a joy to use. In order to make full use of PyBioMed, you will want to know how to write basic programs in Python. Among the many guides to Python, we recommend the documentation at http://www.python.org.


.. _`PyMolecule`: https://raw.githubusercontent.com/gadsbyfly/PyBioMed/master/PyBioMed/download/PyBioMed%20Chem.pdf
.. _`PyProtein`: https://github.com/gadsbyfly/PyBioMed/blob/master/PyBioMed/download/PyBioMed%20Protein.pdf
.. _`PyDNA`: https://github.com/gadsbyfly/PyBioMed/blob/master/PyBioMed/download/PyBioMed%20DNA.pdf
.. _`PyInter`: https://github.com/gadsbyfly/PyBioMed/blob/master/PyBioMed/download/PyBioMed%20Interaction.pdf
