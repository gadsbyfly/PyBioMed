# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
This module is to get different formats of molecules from file and web. If you

have any question please contact me via email.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com
"""

import urllib
import re
import string
import os
from rdkit import Chem

Version=1.0

def ReadMolFromSDF(filename=""):
    """
    Read a set of molecules by SDF file format.
    
    Note: the output of this function is a set of molecular objects.
    
    You need to use for statement to call each object.
    
    Usage:
        
        res=ReadMolFromSDF(filename)
        
        Input: filename is a file name with path.
        
        Output: res is a set of molecular object.
        
    """
    molset=Chem.SDMolSupplier(filename)
    return molset
    


def ReadMolFromMOL(filename=""):
    """
    Read a  molecule by mol file format.
    
    Usage:
        
        res=ReadMolFromMOL(filename)
        
        Input: filename is a file name with path.
        
        Output: res is a  molecular object.
        
    """
    mol=Chem.MolFromMolFile(filename)
    return mol



def ReadMolFromSmile(smi=""):
    """
    #################################################################
    Read a molecule by SMILES string.
        
    Usage:
            
        res=ReadMolFromSmile(smi)
            
        Input: smi is a SMILES string.
            
        Output: res is a molecule object.
    #################################################################
    """
    mol = Chem.MolFromSmiles(string.strip(smi))
        
    return mol
        
        
def ReadMolFromInchi(inchi=""):
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
    mol = Chem.MolFromSmiles(string.strip(smi))
        
    return mol
 
       
def ReadMolFromMol(filename=""):
    """
    #################################################################
    Read a molecule with mol file format.
        
    Usage:
            
        res=ReadMolFromMol(filename)
            
        Input: filename is a file name.
            
        Output: res is a molecule object.
    #################################################################
    """
    mol=Chem.MolFromMolFile(filename)
    return mol
#############################################################################

def GetMolFromCAS(casid=""):
    """
    Downloading the molecules from http://www.chemnet.com/cas/ by CAS ID (casid).
    if you want to use this function, you must be install pybel.
    """
    import pybel
    casid=string.strip(casid)
    localfile=urllib.urlopen('http://www.chemnet.com/cas/supplier.cgi?terms='+casid+'&l=&exact=dict')
    temp=localfile.readlines()
    for i in temp:
        if re.findall('InChI=',i)==['InChI=']:
            k=i.split('    <td align="left">')
            kk=k[1].split('</td>\r\n')
            if kk[0][0:5]=="InChI":
                res=kk[0]    
            else:
                res="None"
    localfile.close()
    mol=pybel.readstring('inchi',string.strip(res))
    smile=mol.write('smi')
    return string.strip(smile)


def GetMolFromEBI():
    """
    """
    pass


def GetMolFromNCBI(cid=""):
    """
    Downloading the molecules from http://pubchem.ncbi.nlm.nih.gov/ by cid (cid).
    """
    cid=string.strip(cid)
    localfile=urllib.urlopen('http://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?cid='+cid+'&disopt=SaveSDF')
    temp=localfile.readlines()   
    f=file("temp.sdf",'w')
    f.writelines(temp)
    f.close()
    localfile.close()
    m=Chem.MolFromMolFile("temp.sdf")
    os.remove("temp.sdf")
    temp=Chem.MolToSmiles(m,isomericSmiles=True)
    return temp


def GetMolFromDrugbank(dbid=""):
    """
    Downloading the molecules from http://www.drugbank.ca/ by dbid (dbid).
    """
    dbid=string.strip(dbid)
    
    localfile=urllib.urlopen('http://www.drugbank.ca/drugs/'+dbid+'.sdf')
    temp=localfile.readlines()   
    f=file("temp.sdf",'w')
    f.writelines(temp)
    f.close()
    localfile.close()
    m=Chem.MolFromMolFile("temp.sdf")
    os.remove("temp.sdf")
    temp=Chem.MolToSmiles(m,isomericSmiles=True)
    return temp



def GetMolFromKegg(kid=""):
    """
    Downloading the molecules from http://www.genome.jp/ by kegg id (kid).
    """
    ID=str(kid)
    localfile=urllib.urlopen('http://www.genome.jp/dbget-bin/www_bget?-f+m+drug+'+ID)
    temp=localfile.readlines() 
    f=file("temp.mol",'w')
    f.writelines(temp)
    f.close()
    localfile.close()
    m=Chem.MolFromMolFile("temp.mol")
    os.remove("temp.mol")
    temp=Chem.MolToSmiles(m,isomericSmiles=True)
    return temp
#############################################################################

if __name__=="__main__":
    print '-'*10+'START'+'-'*10
    print 'Only PyBioMed is successfully installed the code below can be run！'
    from PyBioMed.PyGetMol.GetProtein import timelimited
    @timelimited(10)
    def run_GetMolFromCAS():
        temp=GetMolFromCAS(casid="50-12-4")
        print temp

    @timelimited(10)
    def run_GetMolFromNCBI():
        temp=GetMolFromNCBI(cid="2244")
        print temp

    @timelimited(10)
    def run_GetMolFromDrugbank():
        temp=GetMolFromDrugbank(dbid="DB00133")
        print temp

    @timelimited(10)
    def run_GetMolFromKegg():
        temp=GetMolFromKegg(kid="D02176")
        print temp

    run_GetMolFromCAS()
    print '-'*25
    run_GetMolFromNCBI()
    print '-'*25
    run_GetMolFromDrugbank()
    print '-'*25
    run_GetMolFromKegg()
    print '-'*10+'END'+'-'*10

