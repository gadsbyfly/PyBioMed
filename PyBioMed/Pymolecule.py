# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################

A class used for computing different types of drug descriptors! 

You can freely use and distribute it. If you have any problem, 

you could contact with us timely.

Authors: Dongsheng Cao and Yizeng Liang.

Date: 2012.09.24

Email: oriental-cds@163.com

##############################################################################
"""

from PyBioMed.PyGetMol import Getmol as getmol
from PyMolecule import kappa
from PyMolecule import charge
from PyMolecule import connectivity 
from PyMolecule import constitution 
from PyMolecule import estate
from PyMolecule import geary
from PyMolecule import moe
from PyMolecule import molproperty 
from PyMolecule import moran
from PyMolecule import moreaubroto 
from PyMolecule import topology
from PyMolecule import fingerprint
from PyMolecule import basak
from PyMolecule import bcut
from PyMolecule import cats2d
from PyMolecule import ghosecrippen
from PyMolecule import AtomTypes

from rdkit import Chem
import string

Version=1.0
FingerprintName=['FP2',
                 'FP3',
                 'FP4',
                 'topological',
                 'Estate',
                 'atompairs',
                 'torsions',
                 'morgan',
                 'ECFP2',
                 'ECFP4',
                 'ECFP6',
                 'MACCS',
                 'FCFP2',
                 'FCFP4',
                 'FCFP6',
                 'Pharm2D2point',
                 'Pharm2D3point',
                 'GhoseCrippen',
                 'PubChem']
############################################################################## 

class PyMolecule:
    

    """
    #################################################################
    A PyDrug class used for computing drug descriptors.
    #################################################################
    """
    def __init__(self):
        """
        #################################################################
        constructor of PyMolecule.
        #################################################################
        """
        pass
    
    
    def ReadMolFromMOL(self,filename=""):
        """
        #################################################################
        Read a molecule by SDF or MOL file format.
        
        Usage:
            
            res=ReadMolFromFile(filename)
            
            Input: filename is a file name.
            
            Output: res is a molecule object.
        #################################################################
        """
        self.mol=Chem.MolFromMolMOL(filename)
        return self.mol
    
    
    def ReadMolFromSmile(self,smi=""):
        """
        #################################################################
        Read a molecule by SMILES string.
        
        Usage:
            
            res=ReadMolFromSmile(smi)
            
            Input: smi is a SMILES string.
            
            Output: res is a molecule object.
        #################################################################
        """
        self.mol = Chem.MolFromSmiles(string.strip(smi))
        
        return self.mol
        
        
    def ReadMolFromInchi(self,inchi=""):
        """
        #################################################################
        Read a molecule by Inchi string.
        
        Usage:
            
            res=ReadMolFromInchi(inchi)
            
            Input: inchi is a InChi string.
            
            Output: res is a molecule object.
        #################################################################
        """
        import pybel
        temp=pybel.readstring("inchi",inchi)
        smi=temp.write("smi")
        self.mol = Chem.MolFromSmiles(string.strip(smi))
        
        return self.mol
 
       
    def ReadMolFromMol(self,filename=""):
        """
        #################################################################
        Read a molecule with mol file format.
        
        Usage:
            
            res=ReadMolFromMol(filename)
            
            Input: filename is a file name.
            
            Output: res is a molecule object.
        #################################################################
        """
        self.mol=Chem.MolFromMolFile(filename)
        return self.mol
 
  
   
    def GetMolFromNCBI(self,ID=""):
        """
        #################################################################
        Get a molecule by NCBI id (e.g., 2244).
        
        Usage:
            
            res=GetMolFromNCBI(ID)
            
            Input: ID is a compound ID (CID) in NCBI.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromNCBI(cid=ID)
        return res
 
 
   
    def GetMolFromEBI(self,ID=""):
        """
        #################################################################
        Get a molecule by EBI id.

        Usage:
            
            res=GetMolFromEBI(ID)
            
            Input: ID is a compound identifier in EBI.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromEBI(ID)
        return res
 
 
   
    def GetMolFromCAS(self,ID=""):
        """
        #################################################################
        Get a molecule by kegg id (e.g., 50-29-3).
        
        Usage:
            
            res=GetMolFromCAS(ID)
            
            Input: ID is a CAS identifier.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromCAS(casid=ID)
        return res

   
     
    def GetMolFromKegg(self,ID=""):
        """
        #################################################################
        Get a molecule by kegg id (e.g., D02176).
        
        Usage:
            
            res=GetMolFromKegg(ID)
            
            Input: ID is a compound identifier in KEGG.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromKegg(kid=ID)
        return res

 
 
    def GetMolFromDrugbank(self,ID=""):
        """
        #################################################################
        Get a molecule by drugbank id (e.g.,DB00133).
        
        Usage:
            
            res=GetMolFromDrugbank(ID)
            
            Input: ID is a compound identifier in Drugbank.
            
            Output: res is a SMILES string.
        #################################################################
        """
        res=getmol.GetMolFromDrugbank(dbid=ID)
        return res
  
    
    def GetKappa(self):
        """
        #################################################################
        Calculate all kappa descriptors (7).
        
        Usage:
            
            res=GetKappa()
            
            res is a dict form.
        #################################################################
        """
        res=kappa.GetKappa(self.mol)
        return res
    
    
    def GetCharge(self):
        """
        #################################################################
        Calculate all charge descriptors (25).
        
        Usage:
            
            res=GetCharge()
            
            res is a dict form.
        #################################################################
        """
        res=charge.GetCharge(self.mol)
        return res
    
    
    def GetConnectivity(self):
        """
        #################################################################
        Calculate all conenctivity descriptors (44).
        
        Usage:
            
            res=GetConnectivity()
            
            res is a dict form.
        #################################################################
        """
        res=connectivity.GetConnectivity(self.mol)
        return res
        
    
    
    def GetConstitution(self):
        """
        #################################################################
        Calculate all constitutional descriptors (30).
        
        Usage:
            
            res=GetConstitution()
            
            res is a dict form.
        #################################################################
        """
        res=constitution.GetConstitutional(self.mol)
        return res
        
        
    def GetBasak(self):
        """
        #################################################################
        Calculate all basak's information content descriptors (21).
        
        Usage:
            
            res=GetBasak()
            
            res is a dict form.
        #################################################################
        """
        res=basak.Getbasak(self.mol)
        return res    


    def GetBurden(self):
        
        """
        #################################################################
        Calculate all Burden descriptors (64).
        
        Usage:
            
            res=GetBurden()
            
            res is a dict form.
        #################################################################
        """
        res=bcut.GetBurden(self.mol)
        return res        
        
    
    def GetEstate(self):
        """
        #################################################################
        Calculate estate descriptors (316).
        
        Usage:
            
            res=GetEstate()
            
            res is a dict form.
        #################################################################
        """
        res=estate._GetEstate(self.mol)
        return res
    
    
    def GetGeary(self):
        """
        #################################################################
        Calculate all Geary autocorrelation descriptors (32).
        
        Usage:
            
            res=GetGeary()
            
            res is a dict form.
        #################################################################
        """
        res=geary.GetGearyAuto(self.mol)
        return res
    
    
    def GetMOE(self):
        """
        #################################################################
        Calculate all MOE-type descriptors (60).
        
        Usage:
            
            res=GetMOE()
            
            res is a dict form.
        #################################################################
        """
        res=moe.GetMOE(self.mol)
        return res
    
    
    def GetMolProperty(self):
        """
        #################################################################
        Calculate all molecular properties (6).
        
        Usage:
            
            res=GetMolProperty()
            
            res is a dict form.
        #################################################################
        """
        res=molproperty.GetMolecularProperty(self.mol)
        return res
    
    
    def GetMoran(self):
        """
        #################################################################
        Calculate all Moran autocorrealtion descriptors (32).
        
        Usage:
            
            res=GetMoran()
            
            res is a dict form.
        #################################################################
        """
        res=moran.GetMoranAuto(self.mol)
        return res
    
    
    def GetMoreauBroto(self):
        """
        #################################################################
        Calculate all Moreau-Broto autocorrelation descriptors(32).
        
        Usage:
            
            res=GetMoreauBroto()
            
            res is a dict form.
        #################################################################
        """
        res=moreaubroto.GetMoreauBrotoAuto(self.mol)
        return res
    
    
    def GetTopology(self):
        """
        #################################################################
        Calculate all topological descriptors (25).
        
        Usage:
            
            res=GetTopology()
            
            res is a dict form.
        #################################################################
        """
        res=topology.GetTopology(self.mol)
        return res
    
    
    def GetFingerprint(self,FPName='topological', **kwargs):
        """
        #################################################################
        Calculate all fingerprint descriptors.
        
        see the fingerprint type in FingerprintName
        
        Usage:
            
            res=GetFingerprint(FPName='topological')
            
            res is a tuple or list or dict.
        #################################################################
        """
        
        if FPName in FingerprintName:
            temp=fingerprint._FingerprintFuncs[FPName]
            res=temp(self.mol, **kwargs)
            return res
        else:
            # res=fingerprint.CalculateDaylightFingerprint(self.mol)
            res='This is not a valid fingerprint name！!'
            return res

    
    
    def GetCATS2D(self):
        """
        #################################################################
        The main program for calculating the CATS descriptors.
        
        CATS: chemically advanced template serach
        
        ----> CATS_DA0 ....
        
        Usage:
            
            result=CATS2D(mol,PathLength = 10,scale = 1)
            
            Input: mol is a molecule object.
            
                   PathLength is the max topological distance between two atoms.
                   
                   scale is the normalization method (descriptor scaling method)
                   
                   scale = 1 indicates that no normalization. That is to say: the 
                   
                   values of the vector represent raw counts ("counts").
                   
                   scale = 2 indicates that division by the number of non-hydrogen
                   
                   atoms (heavy atoms) in the molecule.
                   
                   scale = 3 indicates that division of each of 15 possible PPP pairs
                   
                   by the added occurrences of the two respective PPPs.
            
            Output: result is a dict format with the definitions of each descritor.
        #################################################################
        """
        res = cats2d.CATS2D(self.mol, PathLength = 10, scale = 3)        
        
        return res

    # def GetGhoseCrippenFingerprint(self, FPName='GhoseCrippenFingerprint'):
    #     """
    #     #################################################################
    #     Ghose-Crippen substructures based on the definitions of
    #
    #     SMARTS from Ghose-Crippen's paper. (110 dimension)
    #
    #     The result is a dict format.
    #     #################################################################
    #     """
    #     res = ghosecrippen.GhoseCrippenFingerprint(self.mol)
    #
    #     return res
    #
    #
    # def GetGhoseCrippen(self, FPName='GetGhoseCrippen'):
    #     """
    #     #################################################################
    #     Ghose-Crippen counts based on the definitions of
    #
    #     SMARTS from Ghose-Crippen's paper. (110 dimension)
    #
    #     The result is a dict format.
    #     #################################################################
    #     """
    #     res = ghosecrippen.GhoseCrippenFingerprint(self.mol, count = True)
    #
    #     return res

    
    
    
    def GetAllDescriptor(self):
        """
        #################################################################
        Calculate all descriptors (608).
        
        Usage:
            
            res=GetAllDescriptor()
            
            res is a dict form.
        #################################################################
        """
        res={}
        res.update(self.GetKappa())
        res.update(self.GetCharge())
        res.update(self.GetConnectivity())
        res.update(self.GetConstitution())
        res.update(self.GetEstate())
        res.update(self.GetGeary())
        res.update(self.GetMOE())
        res.update(self.GetMoran())
        res.update(self.GetMoreauBroto())
        res.update(self.GetTopology())
        res.update(self.GetMolProperty())
        res.update(self.GetBasak())
        res.update(self.GetBurden())
        res.update(self.GetCATS2D())
        
        return res
##############################################################################    

if __name__=="__main__":
    
    drugclass=PyMolecule()
    drugclass.ReadMolFromSmile("CCC1(c2ccccc2)C(=O)N(C)C(=N1)O")
    print drugclass.GetCharge()    
    
    print drugclass.GetKappa()
    print len(drugclass.GetKappa())
    print drugclass.GetTopology()
    print len(drugclass.GetTopology())
    print drugclass.GetMoreauBroto()
    res=drugclass.GetAllDescriptor()
    print len(res)
    #print drugclass.GetMolFromDrugbank(ID="DB00133")
    #res=drugclass.GetFingerprint(FPName='Estate')
    print res
    print len(res)
    print drugclass.GetConnectivity()
    DrugBankID = 'DB01014'
    drugclass=PyMolecule()
    smi = drugclass.GetMolFromDrugbank(DrugBankID)
    drugclass.ReadMolFromSmile(smi)
    print drugclass.GetKappa()
    
    print drugclass.GetCATS2D()
    print drugclass.GetFingerprint(FPName='Estate')
    # print drugclass.GetGhoseCrippen()
    # print drugclass.GetGhoseCrippenFingerprint()
    print len(drugclass.GetBasak())
    print len(drugclass.GetBurden())
    








   
   
   
   
   
   
   
   
   
   
   
    
    
    
    
    
    