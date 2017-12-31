# -*- coding: utf-8 -*-
"""
##############################################################################

The calculation of Moran autocorrelation descriptors. You can get 32 molecular

decriptors. You can freely use and distribute it. If you hava  any problem, 

you could contact with us timely!

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

def _CalculateMoranAutocorrelation(mol,lag=1,propertylabel='m'):
    """
    #################################################################
    **Internal used only**
    
    Calculation of Moran autocorrelation descriptors based on 
    
    different property weights.
    
    Usage:
        
    res=_CalculateMoranAutocorrelation(mol,lag=1,propertylabel='m')
    
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
                res=res+(temp1-aveweight)*(temp2-aveweight)
                index=index+1
            else:
                res=res+0.0
                
                
    if sum(tempp)==0 or index==0:
        result=0
    else:
        result=(res/index)/(sum(tempp)/Natom)
                
    return round(result,3)


def CalculateMoranAutoMass(mol):
    """
    #################################################################
    Calculation of Moran autocorrelation descriptors based on 
    
    carbon-scaled atomic mass.
    
    Usage:
    
    res=CalculateMoranAutoMass(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight moran autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['MATSm'+str(i+1)]=_CalculateMoranAutocorrelation(mol,lag=i+1,propertylabel='m')
    
    
    return res


def CalculateMoranAutoVolume(mol):
    """
    #################################################################
    Calculation of Moran autocorrelation descriptors based on 
    
    carbon-scaled atomic van der Waals volume.

    Usage:
    
    res=CalculateMoranAutoVolume(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight moran autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['MATSv'+str(i+1)]=_CalculateMoranAutocorrelation(mol,lag=i+1,propertylabel='V')
    
    
    return res

def CalculateMoranAutoElectronegativity(mol):
    """
    #################################################################
    Calculation of Moran autocorrelation descriptors based on 
    
    carbon-scaled atomic Sanderson electronegativity.
    
    Usage:
    
    res=CalculateMoranAutoElectronegativity(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight moran autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['MATSe'+str(i+1)]=_CalculateMoranAutocorrelation(mol,lag=i+1,propertylabel='En')
    
    
    return res

def CalculateMoranAutoPolarizability(mol):
    """
    #################################################################
    Calculation of Moran autocorrelation descriptors based on 
    
    carbon-scaled atomic polarizability.
    
    Usage:
    
    res=CalculateMoranAutoPolarizability(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing eight moran autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    
    for i in range(8):
        res['MATSp'+str(i+1)]=_CalculateMoranAutocorrelation(mol,lag=i+1,propertylabel='alapha')
    
    
    return res


def GetMoranAuto(mol):
    """
    #################################################################
    Calcualate all Moran autocorrelation descriptors.
    
    (carbon-scaled atomic mass, carbon-scaled atomic van der Waals volume,
     
    carbon-scaled atomic Sanderson electronegativity,
     
    carbon-scaled atomic polarizability)
    
    Usage:
    
    res=GetMoranAuto(mol)
    
    Input: mol is a molecule object.
    
    Output: res is a dict form containing all moran autocorrealtion
    
    descriptors.
    #################################################################
    """
    res={}
    res.update(CalculateMoranAutoMass(mol))
    res.update(CalculateMoranAutoVolume(mol))
    res.update(CalculateMoranAutoElectronegativity(mol))
    res.update(CalculateMoranAutoPolarizability(mol))
    
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
        
        print GetMoranAuto(m)
