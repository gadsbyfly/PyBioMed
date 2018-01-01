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

from PyBioMed.PyPretreat.PyPretreatPro import ProteinCheck
from PyBioMed.PyPretreat.PyPretreatDNA import *
from PyBioMed.PyPretreat.PyPretreatMol import *
from rdkit import Chem


def test_pypretreat():


    #==============================================================================
    # PyPretreatDNA
    #==============================================================================
    print '...............................................................'
    print 'testing the PyPretreatDNA module'
        
    
    phyche_index = [[1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455,-0.070, -0.918, 1.603, -0.654, 0.567, 1.019]]
    print (NormalizeIndex(phyche_index,is_convert_dict = False)[0])

    #==============================================================================
    # PyPretreatPro
    #==============================================================================
    print '...............................................................'
    print 'testing the PyPretreatPro module'
    
    protein="ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASU"
    protein="ADGC"
    print ProteinCheck(protein)
    
    #==============================================================================
    # PyPretreatMol
    #==============================================================================
    print '...............................................................'
    print 'testing the PyPretreatMol module'
    
    smiles = ['O=C([O-])c1ccccc1','C[n+]1c([N-](C))cccc1','[2H]C(Cl)(Cl)Cl']
    mol = Chem.MolFromSmiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
    sm = StandardizeMol()
    mol = sm.addhs(mol)
    mol = sm.disconnect_metals(mol)
    mol = sm.largest_fragment(mol)
    mol = sm.normalize(mol)
    mol = sm.uncharge(mol)
    mol = sm.canonicalize_tautomer(mol)
    mol = sm.reionize(mol)
    mol = sm.rmhs(mol)
    mol = sm.addhs(mol)
    print Chem.MolToSmiles(mol, isomericSmiles=True)
    
    #==============================================================================
    # 
    #==============================================================================

if __name__ == '__main__':
    test_pypretreat()


