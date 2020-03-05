# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
This module is to compute the various fingerprints  based on the provided 

fingerprint system. If you have any question please contact me via email.

2016.11.15

@author: Zhijiang Yao and Dongsheng Cao

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import MACCSkeys
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit import DataStructs
from estate import CalculateEstateFingerprint as EstateFingerprint
import pybel
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.Chem.Pharm2D import Generate
from ghosecrippen import GhoseCrippenFingerprint
from PubChemFingerprints import calcPubChemFingerAll

Version=1.0
similaritymeasure=[i[0] for i in DataStructs.similarityFunctions]

################################################################
def CalculateFP2Fingerprint(mol):
    """
    #################################################################
    Calculate FP2 fingerprints (1024 bits).
    
    Usage:
        
        result=CalculateFP2Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger = 1024
    temp = mol.calcfp().bits
    for i in temp:
        res.update({i:1})
    
    return NumFinger,res

def CalculateFP3Fingerprint(mol):
    """
    #################################################################
    Calculate FP3 fingerprints (210 bits).
    
    Usage:
        
        result=CalculateFP3Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger = 210
    temp = mol.calcfp('FP3').bits
    for i in temp:
        res.update({i:1})
    
    return NumFinger,res       

def CalculateFP4Fingerprint(mol):
    """
    #################################################################
    Calculate FP4 fingerprints (307 bits).
    
    Usage:
        
        result=CalculateFP4Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger=307
    temp=mol.calcfp('FP4').bits
    for i in temp:
        res.update({i:1})

    return NumFinger,res
    

def CalculateDaylightFingerprint(mol):
    """
    #################################################################
    Calculate Daylight-like fingerprint or topological fingerprint
    
    (2048 bits).
    
    Usage:
        
        result=CalculateDaylightFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger=2048
    bv=FingerprintMols.FingerprintMol(mol)
    temp=tuple(bv.GetOnBits())
    for i in temp:
        res.update({i:1})

    return NumFinger,res,bv



def CalculateMACCSFingerprint(mol):
    """
    #################################################################
    Calculate MACCS keys (166 bits).
    
    Usage:
        
        result=CalculateMACCSFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res={}
    NumFinger=166
    bv=MACCSkeys.GenMACCSKeys(mol)
    temp=tuple(bv.GetOnBits())
    for i in temp:
        res.update({i:1})

    return NumFinger,res,bv




def CalculateEstateFingerprint(mol):
    """
    #################################################################
    Calculate E-state fingerprints (79 bits).
    
    Usage:
        
        result=CalculateEstateFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    NumFinger=79
    res={}
    temp=EstateFingerprint(mol)
    for i in temp:
        if temp[i]>0:
            res[i[7:]]=1
       
    return NumFinger,res,temp


def CalculateAtomPairsFingerprint(mol):
    """
    #################################################################
    Calculate atom pairs fingerprints
    
    Usage:
        
        result=CalculateAtomPairsFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=Pairs.GetAtomPairFingerprint(mol)
    
    return res.GetLength(),res.GetNonzeroElements(),res



def CalculateTopologicalTorsionFingerprint(mol):
    """
    #################################################################
    Calculate Topological Torsion Fingerprints
    
    Usage:
        
        result=CalculateTopologicalTorsionFingerprint(mol)
        
        Input: mol is a molecule object.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=Torsions.GetTopologicalTorsionFingerprint(mol)

    return res.GetLength(),res.GetNonzeroElements(),res



def CalculateMorganFingerprint(mol,radius=2):
    """
    #################################################################
    Calculate Morgan
    
    Usage:
        
        result=CalculateMorganFingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the number of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=AllChem.GetMorganFingerprint(mol,radius)
    
    return res.GetLength(),res.GetNonzeroElements(),res

def CalculateECFP2Fingerprint(mol,radius=1):
    """
    #################################################################
    Calculate ECFP2
    
    Usage:
        
        result=CalculateECFP2Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the vector of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=AllChem.GetMorganFingerprint(mol,radius)
    
    fp = tuple(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits = 1024))
    
    return fp, res.GetNonzeroElements(), res
    
def CalculateECFP4Fingerprint(mol,radius=2):
    """
    #################################################################
    Calculate ECFP4
    
    Usage:
        
        result=CalculateECFP4Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the vector of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=AllChem.GetMorganFingerprint(mol,radius)
    
    fp = tuple(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits = 1024))
    
    return fp, res.GetNonzeroElements(), res

def CalculateECFP6Fingerprint(mol,radius=3):
    """
    #################################################################
    Calculate ECFP6
    
    Usage:
        
        result=CalculateECFP6Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the vector of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res=AllChem.GetMorganFingerprint(mol,radius)   
    
    fp = tuple(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits = 1024))
    
    return fp, res.GetNonzeroElements(), res
        
def CalculateSimilarityPybel(fp1,fp2):
    """
    #################################################################
    Calculate Tanimoto similarity between two molecules.
    
    Usage:
        
        result=CalculateSimilarityPybel(fp1,fp2)
        
        Input: fp1 and fp2 are two DataStructs.
        
        Output: result is a Tanimoto similarity value.
    #################################################################
    """
    intersection = set(fp1[1].keys())& set(fp2[1].keys())
    union = set(fp1[1].keys()) | set(fp2[1].keys())
    tanimoto = len(intersection) / float(len(union))
    return round(tanimoto,3)


def CalculateSimilarityRdkit(fp1,fp2,similarity="Tanimoto"):
    """
    #################################################################
    Calculate similarity between two molecules.
    
    Usage:
        
        result=CalculateSimilarity(fp1,fp2)
        Users can choose 11 different types:
        Tanimoto, Dice, Cosine, Sokal, Russel,
        RogotGoldberg, AllBit, Kulczynski, 
        McConnaughey, Asymmetric, BraunBlanquet       
        Input: fp1 and fp2 are two DataStructs.
        
        Output: result is a similarity value.
    #################################################################
    """
    temp=DataStructs.similarityFunctions
    for i in temp:
        if similarity in i[0]:
            similarityfunction=i[1]
        else:
            similarityfunction=temp[0][1]
            
    res=similarityfunction(fp1,fp2)
    return round(res,3)


def CalculateFCFP2Fingerprint(mol, radius=1, nBits = 1024):
    """
    #################################################################
    Calculate FCFP2
    
    Usage:
        
        result=CalculateFCFP2Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the vector of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res = AllChem.GetMorganFingerprint(mol, radius, useFeatures = True)
    
    fp = tuple(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits, useFeatures = True))
    
    return fp, res.GetNonzeroElements(), res


def CalculateFCFP4Fingerprint(mol,radius=2, nBits = 1024):
    """
    #################################################################
    Calculate FCFP4
    
    Usage:
        
        result=CalculateFCFP4Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the vector of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res = AllChem.GetMorganFingerprint(mol, radius, useFeatures = True)
    
    fp = tuple(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits, useFeatures = True))
    
    return fp, res.GetNonzeroElements(), res


def CalculateFCFP6Fingerprint(mol,radius=3, nBits = 1024):
    """
    #################################################################
    Calculate FCFP6
    
    Usage:
        
        result=CalculateFCFP4Fingerprint(mol)
        
        Input: mol is a molecule object.
        
        radius is a radius.
        
        Output: result is a tuple form. The first is the vector of 
        
        fingerprints. The second is a dict form whose keys are the 
        
        position which this molecule has some substructure. The third
        
        is the DataStructs which is used for calculating the similarity.
    #################################################################
    """
    res = AllChem.GetMorganFingerprint(mol, radius, useFeatures = True)
    
    fp = tuple(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits, useFeatures = True))
    
    return fp, res.GetNonzeroElements(),res

################################################################
fdefstr = '''
AtomType NDonor [N&!H0&v3,N&!H0&+1&v4,n&H1&+0]
AtomType ChalcDonor [O,S;H1;+0]
DefineFeature SingleAtomDonor [{NDonor},{ChalcDonor},!$([D1]-[C;D3]=[O,S,N])]
  Family Donor
  Weights 1
EndFeature

AtomType NAcceptor [$([N&v3;H1,H2]-[!$(*=[O,N,P,S])])]
Atomtype NAcceptor [$([N;v3;H0])]
AtomType NAcceptor [$([n;+0])]
AtomType ChalcAcceptor [$([O,S;H1;v2]-[!$(*=[O,N,P,S])])]
AtomType ChalcAcceptor [O,S;H0;v2]
Atomtype ChalcAcceptor [O,S;-]
Atomtype ChalcAcceptor [o,s;+0]
AtomType HalogenAcceptor [F]
DefineFeature SingleAtomAcceptor [{NAcceptor},{ChalcAcceptor},{HalogenAcceptor}]
  Family Acceptor
  Weights 1
EndFeature

# this one is delightfully easy:
DefineFeature AcidicGroup [C,S](=[O,S,P])-[O;H1,H0&-1]
  Family NegIonizable
  Weights 1.0,1.0,1.0
EndFeature

AtomType CarbonOrArom_NonCarbonyl [$([C,a]);!$([C,a](=O))]
AtomType BasicNH2 [$([N;H2&+0][{CarbonOrArom_NonCarbonyl}])]
AtomType BasicNH1 [$([N;H1&+0]([{CarbonOrArom_NonCarbonyl}])[{CarbonOrArom_NonCarbonyl}])]
AtomType BasicNH0 [$([N;H0&+0]([{CarbonOrArom_NonCarbonyl}])([{CarbonOrArom_NonCarbonyl}])[{CarbonOrArom_NonCarbonyl}])]
AtomType BasicNakedN [N,n;X2;+0]
DefineFeature BasicGroup [{BasicNH2},{BasicNH1},{BasicNH0},{BasicNakedN}]
  Family PosIonizable
  Weights 1.0
EndFeature

# aromatic rings of various sizes:
DefineFeature Arom5 a1aaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0
EndFeature
DefineFeature Arom6 a1aaaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
DefineFeature Arom7 a1aaaaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
DefineFeature Arom8 a1aaaaaaa1
  Family Aromatic
  Weights 1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0
EndFeature
'''
featFactory = ChemicalFeatures.BuildFeatureFactoryFromString(fdefstr)


def CalculatePharm2D2pointFingerprint(mol, featFactory = featFactory):
    """
    Calculate Pharm2D2point Fingerprints
    """
    sigFactory_2point = SigFactory(featFactory,minPointCount=2,maxPointCount=2)
    sigFactory_2point.SetBins([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)])
    sigFactory_2point.Init()
    res = Generate.Gen2DFingerprint(mol,sigFactory_2point)
    
    res_keys = tuple(res.GetOnBits())
    init_list = [0]*135
    for res_key in res_keys:
        init_list[res_key] = 1
        
    BitVect = tuple(init_list)   
    
    return BitVect, res_keys, res
################################################################
def CalculatePharm2D3pointFingerprint(mol, featFactory = featFactory):
    """
    Calculate Pharm2D3point Fingerprints
    """
    sigFactory_3point = SigFactory(featFactory,minPointCount=3,maxPointCount=3)
    sigFactory_3point.SetBins([(0, 2), (2,4), (4,6), (6,10)])
    sigFactory_3point.Init()
    res = Generate.Gen2DFingerprint(mol,sigFactory_3point)
    
    res_keys = tuple(res.GetOnBits())
    init_list = [0]*2135
    for res_key in res_keys:
        init_list[res_key] = 1
        
    BitVect = tuple(init_list)   
    
    return BitVect, res_keys, res

################################################################
def CalculateGhoseCrippenFingerprint(mol, count = False):
    """
    Calculate GhoseCrippen Fingerprints
    """
    res = GhoseCrippenFingerprint(mol, count=count)
    return res

def CalculatePubChemFingerprint(mol):
    """
    Calculate PubChem Fingerprints
    """
    res = calcPubChemFingerAll(mol)
    return res


_FingerprintFuncs={'FP2':CalculateFP2Fingerprint,
                 'FP3':CalculateFP3Fingerprint,
                 'FP4':CalculateFP4Fingerprint,
                 'topological':CalculateDaylightFingerprint,
                 'Estate':CalculateEstateFingerprint,
                 'atompairs':CalculateAtomPairsFingerprint,
                 'torsions':CalculateTopologicalTorsionFingerprint,
                 'morgan':CalculateMorganFingerprint,
                 'ECFP2':CalculateECFP2Fingerprint,
                 'ECFP4':CalculateECFP4Fingerprint,
                 'ECFP6':CalculateECFP6Fingerprint,
                 'MACCS':CalculateMACCSFingerprint,
                 'FCFP2':CalculateFCFP2Fingerprint,
                 'FCFP4':CalculateFCFP4Fingerprint,
                 'FCFP6':CalculateFCFP6Fingerprint,
                 'Pharm2D2point':CalculatePharm2D2pointFingerprint,
                 'Pharm2D3point':CalculatePharm2D3pointFingerprint,
                   'PubChem': CalculatePubChemFingerprint,
                 'GhoseCrippen': CalculateGhoseCrippenFingerprint}
################################################################


if __name__=="__main__":
    
    print('-'*10+'START'+'-'*10)
    
    ms = [Chem.MolFromSmiles('CCOC=N'), Chem.MolFromSmiles('NC1=NC(=CC=N1)N1C=CC2=C1C=C(O)C=C2')]
    m2 = [pybel.readstring("smi",'CCOC=N'),pybel.readstring("smi",'CCO')]
    res1=CalculateECFP4Fingerprint(ms[0])
    print(res1)
    print('-'*25)
    res2=CalculateECFP4Fingerprint(ms[1])
    print(res2)
    print('-'*25)
    mol = pybel.readstring("smi", 'CCOC=N') 
    res3 = CalculateFP3Fingerprint(mol)
    print(res3)
    print('-'*25)
    
    mol = Chem.MolFromSmiles('O=C1NC(=O)NC(=O)C1(C(C)C)CC=C')
    res4 = CalculatePharm2D2pointFingerprint(mol)[0]
    print(res4)
    print('-'*25)
    res5 = CalculatePharm2D3pointFingerprint(mol)[0]
    print(res5)
    print('-'*25)
    res6 = CalculateGhoseCrippenFingerprint(mol)
    print(res6)
    print('-'*25)
    res7 = CalculatePubChemFingerprint(mol)
    print(res7)
    print('-'*10+'END'+'-'*10)