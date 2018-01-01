# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""

This module is to calculate the ghosecrippen descriptor. If you

have any question please contact me via email.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

"""
import string
import os
from rdkit import Chem

Version = 1.0
###########################################################################
def _ReadPatts(fileName):
    
  """ 
  #################################################################
  *Internal Use Only*
  
  parses the pattern list from the data file
  #################################################################
  """
  patts = {}
  order = []
  with open(fileName,'r') as f:
    lines = f.readlines()
  for line in lines:
    if line[0] != '#':
      splitLine = line.split('\t')
      if len(splitLine) >= 4 and splitLine[0] != '':
        sma = splitLine[1]
        if sma != 'SMARTS':
          sma.replace('"','')
          p = Chem.MolFromSmarts(sma)
          if p:
            cha = string.strip(splitLine[0])
            if cha not in order:
              order.append(cha)
            l = patts.get(cha,[])
            l.append((sma,p))
            patts[cha] = l
        else:
          print('Problems parsing smarts: %s'%(sma))
  return order,patts


###########################################################################
def GhoseCrippenFingerprint(mol,count = False):
    """
    #################################################################
    Ghose-Crippen substructures or counts based on the definitions of 
    
    SMARTS from Ghose-Crippen's paper. (110 dimension)
    
    The result is a dict format.
    #################################################################
    """
    order, patts = _ReadPatts(os.path.dirname(os.path.abspath(__file__))+"/Crippen.txt")
    
    GCres = dict()
    for sma in patts:        
        match = mol.GetSubstructMatches(patts[sma][0][1],False,False)
        temp = len([i[0] for i in match])
        GCres.update({sma:temp})
    
    res = {}
    if count == False:
        for i in GCres:
            if GCres[i] > 0:
                res.update({i:1})
            else:
                res.update({i:0})
    else:
        res =  GCres
        
    return res
                
                
###############################################################################
if __name__ =='__main__':
    smif = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    AllDes = []
    for i in smif:        
        mol=Chem.MolFromSmiles(i)
        order, patts = _ReadPatts(os.path.dirname(os.path.abspath(__file__))+"/Crippen.txt")
        temp = GhoseCrippenFingerprint(mol,count = True)
        AllDes.append(temp)
    print AllDes