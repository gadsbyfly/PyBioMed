# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
This module mainly implements the calculation of MOE-type descriptors, which 

include LabuteASA, TPSA, slogPVSA, MRVSA, PEOEVSA, EstateVSA and VSAEstate, 

respectively (60).

If you have any question about these indices please contact me via email.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

from rdkit import Chem
from rdkit.Chem import MolSurf as MOE 
from rdkit.Chem.EState import EState_VSA as EVSA


Version=1.0
################################################################

def CalculateLabuteASA(mol):
    """
    #################################################################
    Calculation of Labute's Approximate Surface Area (ASA from MOE)
    
    Usage:
        
        result=CalculateLabuteASA(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    res={}
    temp=MOE.pyLabuteASA(mol,includeHs=1)
    res['LabuteASA']=round(temp,3)
    return res

def CalculateTPSA(mol):
    """
    #################################################################
    Calculation of topological polar surface area based on fragments.
    
    Implementation based on the Daylight contrib program tpsa.
    
    Usage:
        
        result=CalculateTPSA(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    res={}
    temp=MOE.TPSA(mol)
    res['MTPSA']=round(temp,3)
    return res

def CalculateSLOGPVSA(mol,bins=None):
    """
    #################################################################
    MOE-type descriptors using LogP contributions and surface 
    
    area contributions.
    
    logpBins=[-0.4,-0.2,0,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6]
    
    You can specify your own bins to compute some descriptors.
    
    Usage:
        
        result=CalculateSLOGPVSA(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################   
    """
    temp=MOE.SlogP_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['slogPVSA'+str(i)]=round(j,3)
    return res


def CalculateSMRVSA(mol,bins=None):
    """
    #################################################################
    MOE-type descriptors using MR contributions and surface 
    
    area contributions.
    
    mrBins=[1.29, 1.82, 2.24, 2.45, 2.75, 3.05, 3.63,3.8,4.0]
    
    You can specify your own bins to compute some descriptors.
    
    Usage:
        
        result=CalculateSMRVSA(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    temp=MOE.SMR_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['MRVSA'+str(i)]=round(j,3)
    return res


def CalculatePEOEVSA(mol,bins=None):
    
    """
    #################################################################
    MOE-type descriptors using partial charges and surface 
    
    area contributions.
    
    chgBins=[-.3,-.25,-.20,-.15,-.10,-.05,0,.05,.10,.15,.20,.25,.30]
    
    You can specify your own bins to compute some descriptors
    
    Usage:
        
        result=CalculatePEOEVSA(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    temp=MOE.PEOE_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['PEOEVSA'+str(i)]=round(j,3)
    return res    


def CalculateEstateVSA(mol,bins=None):
    """
    #################################################################
    MOE-type descriptors using Estate indices and surface area 
    
    contributions.
    
    estateBins=[-0.390,0.290,0.717,1.165,1.540,1.807,2.05,4.69,9.17,15.0] 
    
    You can specify your own bins to compute some descriptors
    
    Usage:
        
        result=CalculateEstateVSA(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    temp=EVSA.EState_VSA_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['EstateVSA'+str(i)]=round(j,3)
    return res


def CalculateVSAEstate(mol,bins=None):
    """
    #################################################################
    MOE-type descriptors using Estate indices and surface 
    
    area contributions.
    
    vsaBins=[4.78,5.00,5.410,5.740,6.00,6.07,6.45,7.00,11.0] 
    
    You can specify your own bins to compute some descriptors
    
    Usage:
        
        result=CalculateVSAEstate(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    temp=EVSA.VSA_EState_(mol,bins,force=1)
    res={}
    for i,j in enumerate(temp):
        res['VSAEstate'+str(i)]=round(j,3)
    return res
    
    
    
def GetMOE(mol):
    """
    #################################################################
    The calculation of MOE-type descriptors (ALL).
    
    Usage:
        
        result=GetMOE(mol)
        
        Input: mol is a molecule object
        
        Output: result is a dict form 
    #################################################################
    """
    result={}
    result.update(CalculateLabuteASA(mol))
    result.update(CalculateTPSA(mol))
    result.update(CalculateSLOGPVSA(mol,bins=None))
    result.update(CalculateSMRVSA(mol,bins=None))
    result.update(CalculatePEOEVSA(mol,bins=None))
    result.update(CalculateEstateVSA(mol,bins=None))
    result.update(CalculateVSAEstate(mol,bins=None))
    return result

#########################################################################

if __name__=="__main__":
    
    
    smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print index+1
        print smi      
        print '\t',GetMOE(m)
        print '\t', len(GetMOE(m))
        

