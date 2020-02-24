# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
This file provides functions to convert descriptors list of multiple molecules (dicts) into CSV
If you have any questions, please feel free to contact us.
E-mail: biomed@csu.edu.cn

@File name: PyPreTools
@author: Jie Dong and Zhijiang Yao
"""
import csv
def DictToCSV(MultiDictList, csvOutPath):
    """
    Convert descriptors list of multiple molecules (dicts) into CSV
    :param MultiDictList: a list contains multiple dicts
    :param csvOutPath: path to save CSV file
    :return: csvOutPath
    """
    try:
        desHeader = MultiDictList[0].keys()
        desContent = []
        for i in MultiDictList:
            temp=[]
            for j in desHeader:
                temp.append(i.get(j))
            desContent.append(temp)
        f = file(csvOutPath, 'w')
        writer = csv.writer(f)
        writer.writerow(tuple(desHeader))
        for k in desContent:
            writer.writerow(tuple(k))
        f.close()
        return csvOutPath
    except Exception as e:
        return str(e)

def ListToCSV(MultiList, csvOutPath,Name='Des'):
    """
    Convert descriptors list of multiple molecules (lists) into CSV
    :param MultiList: a list contains multiple lists
    :param csvOutPath: path to save CSV file
    :return: csvOutPath
    """
    try:
        desHeader = []
        for index in range(len(MultiList[0])):
            desHeader.append(str(Name)+str(index+1))
        desContent = []
        for i in MultiList:
            desContent.append(i)
        f = file(csvOutPath, 'w')
        writer = csv.writer(f)
        writer.writerow(tuple(desHeader))
        for k in desContent:
            writer.writerow(tuple(k))
        f.close()
        return csvOutPath
    except Exception as e:
        return str(e)

def TupleToCSV(MultiTupleList, csvOutPath,Name='Des'):
    """
    Convert descriptors list of multiple molecules (tuple) into CSV
    :param MultiTupleList: a list contains multiple lists
    :param csvOutPath: path to save CSV file
    :return: csvOutPath
    """
    try:
        desHeader = []
        for index in range(len(MultiTupleList[0])):
            desHeader.append(str(Name)+str(index+1))
        desContent = []
        for i in MultiTupleList:
            desContent.append(i)
        f = file(csvOutPath, 'w')
        writer = csv.writer(f)
        writer.writerow(tuple(desHeader))
        for k in desContent:
            writer.writerow(k)
        f.close()
        return csvOutPath
    except Exception as e:
        return str(e)


if __name__=="__main__":
    print('Only PyBioMed is successfully installed the code below can be run！')

    #  uncomment below code as an example to use if you have successfully installed PyBioMed.

    print('-'*10+'START'+'-'*10)
    from rdkit import Chem
    from PyBioMed.PyMolecule.charge import GetCharge
    smis = ['CCCC','CCCCC','CCCCCC','CC(N)C(=O)O','CC(N)C(=O)[O-].[Na+]']
    smi5=['CCCCCC','CCC(C)CC','CC(C)CCC','CC(C)C(C)C','CCCCCN','c1ccccc1N']

    des_list2 = []
    from PyBioMed.PyMolecule.fingerprint import CalculatePubChemFingerprint

    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        des_list2.append(CalculatePubChemFingerprint(m))
    print(des_list2)

    print(ListToCSV(des_list2,'reeeee.csv','pubchem'))
    print('-'*25)

    des_list = []
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        des_list.append(GetCharge(m))
    print(des_list)

    print(DictToCSV(des_list,'reeee.csv'))
    print('-'*25)

    print('-'*10+'END'+'-'*10)