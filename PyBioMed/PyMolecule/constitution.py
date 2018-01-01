# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
The calculation of molecular constitutional indices based on its topological

structure. You can get 30 molecular connectivity descriptors. You can freely

use and distribute it. If you hava  any problem, you could contact with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

from rdkit import Chem
#from rdkit.Chem import rdchem
from rdkit.Chem import Lipinski as LPK

#import math

Version=1.0
#############################################################

def CalculateMolWeight(mol):
    
    """
    #################################################################
    Calculation of molecular weight
    
    Note that not including H
    
    ---->Weight  
    
    Usage:
        
        result=CalculateMolWeight(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
        
    #################################################################
    """
    MolWeight=0
    for atom in mol.GetAtoms():
        MolWeight=MolWeight+atom.GetMass()

    return MolWeight


def CalculateAverageMolWeight(mol):
    """
    #################################################################
    Calcualtion of average molecular weight
    
    Note that not including H
    
    ---->AWeight
    
    Usage:
        
        result=CalculateAverageMolWeight(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    
    MolWeight=0
    for atom in mol.GetAtoms():
        MolWeight=MolWeight+atom.GetMass()

    return MolWeight/mol.GetNumAtoms()


def CalculateHydrogenNumber(mol):

    """
    #################################################################
    Calculation of Number of Hydrogen in a molecule
    
    ---->nhyd
    
    Usage:
        
        result=CalculateHydrogenNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0
    Hmol=Chem.AddHs(mol)
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==1:
            i=i+1
            
    return i

def CalculateHalogenNumber(mol):

    """
    #################################################################
    Calculation of Halogen counts in a molecule
    
    ---->nhal
    
    Usage:
        
        result=CalculateHalogenNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==9 or atom.GetAtomicNum()==17 or atom.GetAtomicNum()==35 or atom.GetAtomicNum()==53:
            i=i+1
    return i
            


def CalculateHeteroNumber(mol):
    """
    #################################################################
    Calculation of Hetero counts in a molecule
    
    ---->nhet
    
    Usage:
        
        result=CalculateHeteroNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==6 or atom.GetAtomicNum()==1:
            i=i+1

    return mol.GetNumAtoms()-i



def CalculateHeavyAtomNumber(mol):
    """
    #################################################################
    Calculation of Heavy atom counts in a molecule
    
    ---->nhev
    
    Usage:
        
        result=CalculateHeavyAtomNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return mol.GetNumAtoms(onlyHeavy=1)


def _CalculateElementNumber(mol,AtomicNumber=6):

    """
    #################################################################
    **Internal used only**
    
    Calculation of element counts with atomic number equal to n in a molecule
    #################################################################
    """
    i=0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNumber:
            i=i+1
            
    return i


def CalculateFluorinNumber(mol):

    """
    #################################################################
    Calculation of Fluorin counts in a molecule
    
    ---->ncof
    
    Usage:
        
        result=CalculateFluorinNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
            
    return _CalculateElementNumber(mol,AtomicNumber=9)


def CalculateChlorinNumber(mol):

    """
    #################################################################
    Calculation of Chlorin counts in a molecule
    
    ---->ncocl
    
    Usage:
        
        result=CalculateChlorinNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=17)


def CalculateBromineNumber(mol):

    """
    #################################################################
    Calculation of Bromine counts in a molecule
    
    ---->ncobr
    
    Usage:
        
        result=CalculateBromineNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=35)

def CalculateIodineNumber(mol):
    """
    #################################################################
    Calculation of Iodine counts in a molecule
    
    ---->ncoi
    
    Usage:
        
        result=CalculateIodineNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=53)


def CalculateCarbonNumber(mol):

    """
    #################################################################
    Calculation of Carbon number in a molecule
    
    ---->ncarb
    
    Usage:
        
        result=CalculateCarbonNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """

    return _CalculateElementNumber(mol,AtomicNumber=6)

def CalculatePhosphorNumber(mol):
    """
    #################################################################
    Calcualtion of Phosphor number in a molecule
    
    ---->nphos
    
    Usage:
        
        result=CalculatePhosphorNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementNumber(mol,AtomicNumber=15)


def CalculateSulfurNumber(mol):
    """
    #################################################################
    Calculation of Sulfur counts in a molecule
    
    ---->nsulph
    
    Usage:
        
        result=CalculateSulfurNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementNumber(mol,AtomicNumber=16)



def CalculateOxygenNumber(mol):
    """
    #################################################################
    Calculation of Oxygen counts in a molecule
    
    ---->noxy
    
    Usage:
        
        result=CalculateOxygenNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################

    """
    return _CalculateElementNumber(mol,AtomicNumber=8)
        


def CalculateNitrogenNumber(mol):
    """
    #################################################################
    Calculation of Nitrogen counts in a molecule
    
    ---->nnitro
    
    Usage:
        
        result=CalculateNitrogenNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    
    return _CalculateElementNumber(mol,AtomicNumber=7)


def CalculateRingNumber(mol):
    """
    #################################################################
    Calculation of ring counts in a molecule
    
    ---->nring
    
    Usage:
        
        result=CalculateRingNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return Chem.GetSSSR(mol)


def CalculateRotationBondNumber(mol):
    """
    #################################################################
    Calculation of rotation bonds counts in a molecule
    
    ---->nrot
    
    Note that this is the same as calculation of single bond
    
    counts in a molecule.
    
    Usage:
        
        result=CalculateRotationBondNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return LPK.NumRotatableBonds(mol)


def CalculateHdonorNumber(mol):
    """
    #################################################################
    Calculation of Hydrongen bond donor counts in a molecule
    
    ---->ndonr
    
    Usage:
        
        result=CalculateHdonorNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return LPK.NumHDonors(mol)



def CalculateHacceptorNumber(mol):
    """
    #################################################################
    Calculation of Hydrogen bond acceptor counts in a molecule
    
    ---->naccr
    
    Usage:
        
        result=CalculateHacceptorNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return LPK.NumHAcceptors(mol)


def CalculateSingleBondNumber(mol):

    """
    #################################################################
    Calculation of single bond counts in a molecule
    
    ---->nsb
    
    Usage:
        
        result=CalculateSingleBondNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0;
    for bond in mol.GetBonds():

        if bond.GetBondType().name=='SINGLE':
            i=i+1
            
    return i
            

def CalculateDoubleBondNumber(mol):

    """
    #################################################################
    Calculation of double bond counts in a molecule
    
    ---->ndb
    
    Usage:
        
        result=CalculateDoubleBondNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0;
    for bond in mol.GetBonds():

        if bond.GetBondType().name=='DOUBLE':
            i=i+1
            
    return i
            
        

def CalculateTripleBondNumber(mol):

    """
    #################################################################
    Calculation of triple bond counts in a molecule
    
    ---->ntb
    
    Usage:
        
        result=CalculateTripleBondNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0;
    for bond in mol.GetBonds():

        if bond.GetBondType().name=='TRIPLE':
            i=i+1
            
    return i
            
def CalculateAromaticBondNumber(mol):

    """
    #################################################################
    Calculation of aromatic bond counts in a molecule
    
    ---->naro
    
    Usage:
        
        result=CalculateAromaticBondNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    i=0;
    for bond in mol.GetBonds():

        if bond.GetBondType().name=='AROMATIC':
            i=i+1
            
    return i
    
def CalculateAllAtomNumber(mol):
    """
    #################################################################
    Calculation of all atom counts in a molecule
    
    ---->nta
    
    Usage:
        
        result=CalculateAllAtomNumber(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    
    return Chem.AddHs(mol).GetNumAtoms()
        
def _CalculatePathN(mol,PathLength=2):
    """
    #################################################################
    *Internal Use Only*
    
    Calculation of the counts of path length N for a molecule
    
    ---->PC1-PC6
    
    Usage:
        
        result=CalculateMolWeight(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a numeric value.
    #################################################################
    """
    return len(Chem.FindAllPathsOfLengthN(mol,PathLength,useBonds=1))

def CalculatePath1(mol):
    """
    #################################################################
    Calculation of the counts of path length 1 for a molecule
    #################################################################
    """
    return _CalculatePathN(mol,1)

def CalculatePath2(mol):
    """
    #################################################################
    Calculation of the counts of path length 2 for a molecule
    #################################################################
    """
    return _CalculatePathN(mol,2)

def CalculatePath3(mol):
    """
    #################################################################
    Calculation of the counts of path length 3 for a molecule
    #################################################################
    """
    return _CalculatePathN(mol,3)

def CalculatePath4(mol):
    """
    #################################################################
    Calculation of the counts of path length 4 for a molecule
    #################################################################
    """
    return _CalculatePathN(mol,4)

def CalculatePath5(mol):
    """
    #################################################################
    Calculation of the counts of path length 5 for a molecule
    #################################################################
    """
    return _CalculatePathN(mol,5)

def CalculatePath6(mol):
    """
    #################################################################
    Calculation of the counts of path length 6 for a molecule
    #################################################################
    """
    return _CalculatePathN(mol,6)


_constitutional={'Weight':CalculateMolWeight,
                'AWeight':CalculateAverageMolWeight,
                'nhyd':CalculateHydrogenNumber,
                'nhal':CalculateHalogenNumber,
                'nhet':CalculateHeteroNumber,
                'nhev':CalculateHeavyAtomNumber,
                'ncof':CalculateFluorinNumber,
                'ncocl':CalculateChlorinNumber,
                'ncobr':CalculateBromineNumber,
                'ncoi':CalculateIodineNumber,
                'ncarb':CalculateCarbonNumber,
                'nphos':CalculatePhosphorNumber,
                'nsulph':CalculateOxygenNumber,
                'noxy':CalculateOxygenNumber,
                'nnitro':CalculateNitrogenNumber,
                'nring':CalculateRingNumber,
                'nrot':CalculateRotationBondNumber,
                'ndonr':CalculateHdonorNumber,
                'naccr':CalculateHacceptorNumber,
                'nsb':CalculateSingleBondNumber,
                'ndb':CalculateDoubleBondNumber,
                'naro':CalculateAromaticBondNumber,
                'ntb':CalculateTripleBondNumber,
                'nta':CalculateAllAtomNumber,
                'PC1':CalculatePath1,
                'PC2':CalculatePath2,
                'PC3':CalculatePath3,
                'PC4':CalculatePath4,
                'PC5':CalculatePath5,
                'PC6':CalculatePath6
    }

def GetConstitutional(mol):
    """
    #################################################################
    Get the dictionary of constitutional descriptors for given moelcule mol
    
    Usage:
        
        result=GetConstitutional(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a dict form containing all constitutional values.
    #################################################################
    """
    result={}
    for DesLabel in _constitutional.keys():
        result[DesLabel]=round(_constitutional[DesLabel](mol),3)
    return result

#############################################################

if __name__ =='__main__':
    

    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print index+1
        print smi      
        print '\t',GetConstitutional(m)
        print len(GetConstitutional(m))
    
