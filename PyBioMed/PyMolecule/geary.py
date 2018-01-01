# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
The calculation of Geary autocorrelation indices based on its topological

structure. You can get 32 molecular autocorrelation descriptors. You can 

freely use and distribute it. If you hava  any problem, you could contact 

with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

from rdkit import Chem
from AtomProperty import GetRelativeAtomicProperty

import numpy

Version=1.0


################################################################
def _CalculateGearyAutocorrelation(mol,lag=1,propertylabel='m'):
    """
    #################################################################
    **Internal used only**
    
    Calculation of Geary autocorrelation descriptors based on 
    
    different property weights.
    
    Usage:
        
    res=_CalculateGearyAutocorrelation(mol,lag=1,propertylabel='m')
    
    Input: mol is a molecule object.
    
    lag is the topological distance between atom i and atom j.
    
    propertylabel is the weighted property.
    
    Output: res is a numeric value.
    #################################################################
    """

    Natom=mol.GetNumAtoms()
    
    prolist=[]
    for i in mol.GetAtoms():
        temp=GetRelativeAtomicProperty(i.GetSymbol(),propertyname=propertylabel)
        prolist.append(temp)
        
    aveweight=sum(prolist)/Natom
    
    tempp=[numpy.square(x-aveweight) for x in prolist]   
    
    GetDistanceMatrix=Chem.GetDistanceMatrix(mol)
    res=0.0
    index=0
    for i in range(Natom):
        for j in range(Natom):  
            if GetDistanceMatrix[i,j]==lag:
                atom1=mol.GetAtomWithIdx(i)
                atom2=mol.GetAtomWithIdx(j)
                temp1=GetRelativeAtomicProperty(element=atom1.GetSymbol(),propertyname=propertylabel)
                temp2=GetRelativeAtomicProperty(element=atom2.GetSymbol(),propertyname=propertylabel)
                res=res+numpy.square(temp1-temp2)
                index=index+1
            else:
                res=res+0.0
                
                
    if sum(tempp)==0 or index==0:
        result=0
    else:
        result=(res/index/2)/(sum(tempp)/(Natom-1))
                
    return round(result,3)


def CalculateGearyAutoMass(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic mass.
    
    Usage:
    
    res=CalculateMoranAutoMass(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSm'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='m')
    
    
    return res


def CalculateGearyAutoVolume(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic van der Waals volume.

    Usage:
    
    res=CalculateGearyAutoVolume(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSv'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='V')
    
    
    return res

def CalculateGearyAutoElectronegativity(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic Sanderson electronegativity.
    
    Usage:
    
    res=CalculateGearyAutoElectronegativity(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSe'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='En')
    
    
    return res

def CalculateGearyAutoPolarizability(mol):
    """
    #################################################################
    Calculation of Geary autocorrelation descriptors based on 
    
    carbon-scaled atomic polarizability.
    
    Usage:
    
    res=CalculateGearyAutoPolarizability(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['GATSp'+str(i+1)]=_CalculateGearyAutocorrelation(mol,lag=i+1,propertylabel='alapha')
    
    
    return res


def GetGearyAuto(mol):
    """
    #################################################################
    Calcualate all Geary autocorrelation descriptors.

    (carbon-scaled atomic mass, carbon-scaled atomic van der Waals volume,
     
    carbon-scaled atomic Sanderson electronegativity,
     
    carbon-scaled atomic polarizability)
    
    Usage:
    
    res=GetGearyAuto(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing all geary autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    res.update(CalculateGearyAutoMass(mol))
    res.update(CalculateGearyAutoVolume(mol))
    res.update(CalculateGearyAutoElectronegativity(mol))
    res.update(CalculateGearyAutoPolarizability(mol))
    
    return res
###########################################################################
if __name__=='__main__':
    
    
    smi5=['COCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCOCCN','c1ccccc1N']
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print index+1
        print smi      
##        print '\t',CalculateEstateFingerprint(m)
##        print '\t',CalculateEstateValue(m)
##        print '\t',CalculateMaxAtomTypeEState(m)
##        print '\t', CalculateMinAtomTypeEState(m)
        
        print GetGearyAuto(m)
