# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################

The calculation of Charge descriptors based on Gasteiger/Marseli partial 

charges(25). You can freely use and distribute it. If you hava  any problem, 

you could contact with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

##############################################################################
"""
from rdkit import Chem
from rdkit.Chem import rdPartialCharges as GMCharge

import numpy

Version=1.0
##############################################################################
iter_step=12

def _CalculateElementMaxPCharge(mol,AtomicNum=6):
    """
    #################################################################
    **Internal used only**
    
    Most positive charge on atom with atomic number equal to n
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        return round(max(res),3)

def _CalculateElementMaxNCharge(mol,AtomicNum=6):
    """
    #################################################################
    **Internal used only**
    
    Most negative charge on atom with atomic number equal to n
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(min(res),3)


def CalculateHMaxPCharge(mol):
    """
    #################################################################
    Most positive charge on H atoms
    
    -->QHmax
    
    Usage:
    
        result=CalculateHMaxPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxPCharge(mol,AtomicNum=1)


def CalculateCMaxPCharge(mol):
    """
    #################################################################
    Most positive charge on C atoms
    
    -->QCmax

    Usage:
    
        result=CalculateCMaxPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxPCharge(mol,AtomicNum=6)


def CalculateNMaxPCharge(mol):
    """
    #################################################################
    Most positive charge on N atoms
    
    -->QNmax

    Usage:
    
        result=CalculateNMaxPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxPCharge(mol,AtomicNum=7)


def CalculateOMaxPCharge(mol):
    """
    #################################################################
    Most positive charge on O atoms
    
    -->QOmax

    Usage:
    
        result=CalculateOMaxPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxPCharge(mol,AtomicNum=8)

def CalculateHMaxNCharge(mol):
    """
    #################################################################
    Most negative charge on H atoms
  
    -->QHmin

    Usage:
    
        result=CalculateHMaxNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxNCharge(mol,AtomicNum=1)


def CalculateCMaxNCharge(mol):
    """
    #################################################################
    Most negative charge on C atoms
    
    -->QCmin

    Usage:
    
        result=CalculateCMaxNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxNCharge(mol,AtomicNum=6)


def CalculateNMaxNCharge(mol):
    """
    #################################################################
    Most negative charge on N atoms
    
    -->QNmin

    Usage:
    
        result=CalculateNMaxNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxNCharge(mol,AtomicNum=7)


def CalculateOMaxNCharge(mol):
    """
    #################################################################
    Most negative charge on O atoms
    
    -->QOmin

    Usage:
    
        result=CalculateOMaxNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementMaxNCharge(mol,AtomicNum=8)

def CalculateAllMaxPCharge(mol):
    """
    #################################################################
    Most positive charge on ALL atoms
    
    -->Qmax

    Usage:
    
        result=CalculateAllMaxPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(max(res),3)


def CalculateAllMaxNCharge(mol):
    """
    #################################################################
    Most negative charge on all atoms
    
    -->Qmin

    Usage:
    
        result=CalculateAllMaxNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(min(res),3)


def _CalculateElementSumSquareCharge(mol,AtomicNum=6):
    """
    #################################################################
    **Internal used only**
    
    Ths sum of square Charges on all atoms with atomicnumber equal to n
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
        if atom.GetAtomicNum()==AtomicNum:
            res.append(float(atom.GetProp('_GasteigerCharge')))
    if res==[]:
        return 0
    else:
        return round(sum(numpy.square(res)),3)


def CalculateHSumSquareCharge(mol):
    
    """
    #################################################################
    The sum of square charges on all H atoms
    
    -->QHss

    Usage:
    
        result=CalculateHSumSquareCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementSumSquareCharge(mol,AtomicNum=1)


def CalculateCSumSquareCharge(mol):
    """
    #################################################################
    The sum of square charges on all C atoms
    
    -->QCss

    Usage:
    
        result=CalculateCSumSquareCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementSumSquareCharge(mol,AtomicNum=6)


def CalculateNSumSquareCharge(mol):
    """
    #################################################################
    The sum of square charges on all N atoms
    
    -->QNss

    Usage:
    
        result=CalculateNSumSquareCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementSumSquareCharge(mol,AtomicNum=7)

def CalculateOSumSquareCharge(mol):
    """
    #################################################################
    The sum of square charges on all O atoms
    
    -->QOss

    Usage:
    
        result=CalculateOSumSquareCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    return _CalculateElementSumSquareCharge(mol,AtomicNum=8)

def CalculateAllSumSquareCharge(mol):
    """
    #################################################################
    The sum of square charges on all atoms
    
    -->Qass

    Usage:
    
        result=CalculateAllSumSquareCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        return round(sum(numpy.square(res)),3)

def CalculateTotalPCharge(mol):
    """
    #################################################################
    The total postive charge
    
    -->Tpc

    Usage:
    
        result=CalculateTotalPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(sum(cc[cc>0]),3)

def CalculateMeanPCharge(mol):
    """
    #################################################################
    The average postive charge
    
    -->Mpc
    
    Usage:
    
        result=CalculateMeanPCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(numpy.mean(cc[cc>0]),3)


def CalculateTotalNCharge(mol):
    """
    #################################################################
    The total negative charge
    
    -->Tnc
    
    Usage:
    
        result=CalculateTotalNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(sum(cc[cc<0]),3)


def CalculateMeanNCharge(mol):
    """
    #################################################################
    The average negative charge
    
    -->Mnc
    
    Usage:
    
        result=CalculateMeanNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(numpy.mean(cc[cc<0]),3)


def CalculateTotalAbsoulteCharge(mol):
    """
    #################################################################
    The total absolute charge
    
    -->Tac
    
    Usage:
    
        result=CalculateTotalAbsoulteCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(sum(numpy.absolute(cc)),3)

def CalculateMeanAbsoulteCharge(mol):
    """
    #################################################################
    The average absolute charge
    
    -->Mac
    
    Usage:
    
        result=CalculateMeanAbsoulteCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        return round(numpy.mean(numpy.absolute(cc)),3)

def CalculateRelativePCharge(mol):
    """
    #################################################################
    The partial charge of the most positive atom divided by
    
    the total positive charge.
    
    -->Rpc
    
    Usage:
    
        result=CalculateRelativePCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        if sum(cc[cc>0])==0:
            return 0
        else:
            return round(max(res)/sum(cc[cc>0]),3)

def CalculateRelativeNCharge(mol):
    """
    #################################################################
    The partial charge of the most negative atom divided
    
    by the total negative charge.
    
    -->Rnc
    
    Usage:
    
        result=CalculateRelativeNCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    GMCharge.ComputeGasteigerCharges(Hmol,iter_step)
    res=[]
    for atom in Hmol.GetAtoms():
            res.append(float(atom.GetProp('_GasteigerCharge')))
            
    if res==[]:
        return 0
    else:
        cc=numpy.array(res,'d')
        if sum(cc[cc<0])==0:
            return 0
        else:
            return round(min(res)/sum(cc[cc<0]),3)

def CalculateLocalDipoleIndex(mol):
    """
    #################################################################
    Calculation of local dipole index (D)
    
    -->LDI
    
    Usage:
    
        result=CalculateLocalDipoleIndex(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """

    GMCharge.ComputeGasteigerCharges(mol,iter_step)
    res=[]
    for atom in mol.GetAtoms():
        res.append(float(atom.GetProp('_GasteigerCharge')))
    cc = [numpy.absolute(res[x.GetBeginAtom().GetIdx()]-res[x.GetEndAtom().GetIdx()]) for x in mol.GetBonds()]
    B=len(mol.GetBonds())
    
    return round(sum(cc)/B,3)
        
def CalculateSubmolPolarityPara(mol):
    """
    #################################################################
    Calculation of submolecular polarity parameter(SPP)
    
    -->SPP
    
    Usage:
    
        result=CalculateSubmolPolarityPara(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a numeric value.
    #################################################################
    """

    return round(CalculateAllMaxPCharge(mol)-CalculateAllMaxNCharge(mol),3)


_Charge={'SPP':CalculateSubmolPolarityPara,
        'LDI':CalculateLocalDipoleIndex,
        'Rnc':CalculateRelativeNCharge,
        'Rpc':CalculateRelativePCharge,
        'Mac':CalculateMeanAbsoulteCharge,
        'Tac':CalculateTotalAbsoulteCharge,
        'Mnc':CalculateMeanNCharge,
        'Tnc':CalculateTotalNCharge,
        'Mpc':CalculateMeanPCharge,
        'Tpc':CalculateTotalPCharge,
        'Qass':CalculateAllSumSquareCharge,
        'QOss':CalculateOSumSquareCharge,
        'QNss':CalculateNSumSquareCharge,
        'QCss':CalculateCSumSquareCharge,
        'QHss':CalculateHSumSquareCharge,
        'Qmin':CalculateAllMaxNCharge,
        'Qmax':CalculateAllMaxPCharge,
        'QOmin':CalculateOMaxNCharge,
        'QNmin':CalculateNMaxNCharge,
        'QCmin':CalculateCMaxNCharge,
        'QHmin':CalculateHMaxNCharge,
        'QOmax':CalculateOMaxPCharge,
        'QNmax':CalculateNMaxPCharge,
        'QCmax':CalculateCMaxPCharge,
        'QHmax':CalculateHMaxPCharge,
    }


def GetCharge(mol):
    """
    #################################################################
    Get the dictionary of constitutional descriptors for given moelcule mol
    
    Usage:
    
        result=GetCharge(mol)
    
        Input: mol is a molecule object.
    
        Output: result is a dict form containing all charge descriptors.
    #################################################################
    """
    result={}
    for DesLabel in _Charge.keys():
        result[DesLabel]=_Charge[DesLabel](mol)
    return result

##############################################################################

if __name__ =='__main__':
    

    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index+1)
        print(smi)      
        print('\t',GetCharge(m))
        print(len(GetCharge(m)))


