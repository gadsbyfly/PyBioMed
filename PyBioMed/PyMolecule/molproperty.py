# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
Calculation of Molecular physical/chemical properties based on some special 

type of approaches(6), including: LogP; LogP2; MR; TPSA, UI and Hy.You can 

freely use and distribute it. If you hava  any problem, you could contact 

with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""
from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import MolSurf as MS

import math


Version=1.0
##############################################################
def CalculateMolLogP(mol):
    """
    #################################################################
    Cacluation of LogP value based on Crippen method
    
    ---->LogP
    
    Usage:
        
        result=CalculateMolLogP(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return round(Crippen._pyMolLogP(mol),3)

def CalculateMolLogP2(mol):
    """
    #################################################################
    Cacluation of LogP^2 value based on Crippen method
    
    ---->LogP2
    
    Usage:
        
        result=CalculateMolLogP2(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    res=Crippen._pyMolLogP(mol)
    
    return round(res**2,3)

def CalculateMolMR(mol):
    """
    #################################################################
    Cacluation of molecular refraction value based on Crippen method
    
    ---->MR
    
    Usage:
        
        result=CalculateMolMR(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return round(Crippen._pyMolMR(mol),3)

def CalculateTPSA(mol):
    """
    #################################################################
    calculates the polar surface area of a molecule based upon fragments

    Algorithm in:
        
    P. Ertl, B. Rohde, P. Selzer
    
    Fast Calculation of Molecular Polar Surface Area as a Sum of 
     
    Fragment-based Contributions and Its Application to the Prediction
     
    of Drug Transport Properties, J.Med.Chem. 43, 3714-3717, 2000

    Implementation based on the Daylight contrib program tpsa.
    
    ---->TPSA
    
    Usage:
        
        result=CalculateTPSA(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return round(MS.TPSA(mol),3)


def _CalculateBondNumber(mol,bondtype='SINGLE'):

    """
    ################################################################# 
    **Internal used only*
    
    Calculation of bond counts in a molecule. it may be 
    
    SINGLE, DOUBLE, TRIPLE and AROMATIC
    #################################################################
    """
    i=0;
    for bond in mol.GetBonds():

        if bond.GetBondType().name==bondtype:
            i=i+1
            
    return i


def CalculateUnsaturationIndex(mol):
    """
    #################################################################
    Calculation of unsaturation index.
    
    ---->UI
    
    Usage:
        
        result=CalculateUnsaturationIndex(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    nd=_CalculateBondNumber(mol,bondtype='DOUBLE')
    nt=_CalculateBondNumber(mol,bondtype='TRIPLE')
    na=_CalculateBondNumber(mol,bondtype='AROMATIC')
    res=math.log((1+nd+nt+na),2)
    
    return round(res,3)
    

def CalculateHydrophilicityFactor(mol):
    """
    #################################################################
    Calculation of hydrophilicity factor. The hydrophilicity 
    
    index is described in more detail on page 225 of the 
    
    Handbook of Molecular Descriptors (Todeschini and Consonni 2000).
    
    ---->Hy
    
    Usage:
        
        result=CalculateHydrophilicityFactor(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    nheavy=mol.GetNumAtoms(onlyHeavy=1)
    nc=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==6:
            nc=nc+1
    nhy=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==7 or atom.GetAtomicNum()==8 or atom.GetAtomicNum()==16:
            atomn=atom.GetNeighbors()
            for i in atomn:
                if i.GetAtomicNum()==1:
                    nhy=nhy+1
                
    res=(1+nhy)*math.log((1+nhy),2)+nc*(1.0/nheavy*math.log(1.0/nheavy,2))+math.sqrt((nhy+0.0)/(nheavy**2))
    return round(res,3)
    

def CalculateXlogP(mol):
    """
    #################################################################
    Calculation of Wang octanol water partition coefficient.
    
    ---->XLogP
    
    Usage:
        
        result=CalculateXlogP(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    pass


def CalculateXlogP2(mol):
    """
    #################################################################
    Calculation of Wang octanol water partition coefficient (XLogP^2).
    
    ---->XLogP2
    
    Usage:
        
        result=CalculateMolLogP(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    pass


MolecularProperty={'LogP':CalculateMolLogP,
                   'LogP2':CalculateMolLogP2,
                   'MR':CalculateMolMR,
                   'TPSA':CalculateTPSA,
                   'Hy':CalculateHydrophilicityFactor,
                   'UI':CalculateUnsaturationIndex
    }

def GetMolecularProperty(mol):
    """
    #################################################################
    Get the dictionary of constitutional descriptors for 
    
    given moelcule mol
    
    Usage:
        
        result=GetMolecularProperty(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing 6 molecular properties.
    #################################################################
    """
    result={}
    for DesLabel in MolecularProperty.keys():
        result[DesLabel]=MolecularProperty[DesLabel](mol)
    return result
    
##########################################################

if __name__ =='__main__':

    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index+1)
        print(smi)      
        print('\t',GetMolecularProperty(m))
    #f.close()
