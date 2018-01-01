# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
The script is used for testing.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com 
"""

    
import os

def test_pyinteration():
    
    from PyBioMed.PyInteraction.PyInteraction import CalculateInteraction1
    from PyBioMed.PyInteraction.PyInteraction import CalculateInteraction2
    from PyBioMed.PyInteraction.PyInteraction import CalculateInteraction3
    
        
    
    from PyBioMed.PyDNA import PyDNAac
    
    print '...............................................................'
    print 'testing the DNA descriptors'
    
    
    DNA_des = PyDNAac.GetTCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome','MW-kg'])
    
    print DNA_des
    
    
    
    print '...............................................................'
    print 'testing the protein descriptors'
    
    
    from PyBioMed.PyProtein import CTD
    protein="ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    protein_des = CTD.CalculateCTD(protein)
    
    
    print '...............................................................'
    print 'testing the molecular descriptors'
    
    from PyBioMed.PyMolecule import moe
    from rdkit import Chem
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    m = Chem.MolFromSmiles(smis[3])
    mol_des = moe.GetMOE(m)
    
    
    print '...............................................................'
    print 'testing the Interaction type 1 module'
    
    mol_mol_interaction1 = CalculateInteraction1(mol_des,mol_des)
    print mol_mol_interaction1
    
    pro_mol_interaction1 = CalculateInteraction1(mol_des,protein_des)
    print pro_mol_interaction1
    
    DNA_mol_interaction1 = CalculateInteraction1(DNA_des,mol_des)
    print DNA_mol_interaction1
    
    print '...............................................................'
    print 'testing the Interaction type 2 module'
    
    mol_mol_interaction2 = CalculateInteraction2(mol_des,mol_des)
    print mol_mol_interaction2
    
    pro_mol_interaction2 = CalculateInteraction2(mol_des,protein_des)
    print pro_mol_interaction2
    
    DNA_mol_interaction2 = CalculateInteraction2(DNA_des,mol_des)
    print DNA_mol_interaction2
    
    print '...............................................................'
    print 'testing the Interaction type 3 module'
    
    mol_mol_interaction3 = CalculateInteraction3(mol_des,mol_des)
    print mol_mol_interaction3

if __name__ == '__main__':
    test_pyinteration()





