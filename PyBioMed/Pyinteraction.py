# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
The calculation of interaction descriptors. You can choose three types of

interacation descriptors. You can freely use and distribute it. If you

hava any problem, you could contact with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com
"""

##############################################################################
def CalculateInteraction1(dict1={}, dict2={}):
    """
    Calculate the two interaction features by combining two different

    features.

    Usage:

        res=CalculateInteraction(dict1,dict2)

        Input: dict1 is a dict form containing features.

               dict2 is a dict form containing features.

        Output: res is a dict form containing interaction

        features.
    """
    result = {}
    result.update(dict1)
    for i in dict2:
        result[i + "ex"] = dict2[i]

    return result


def CalculateInteraction2(dict1={}, dict2={}):
    """
    Calculate the two interaction features by combining two different

    features.

    Usage:

        res=CalculateInteraction(dict1,dict2)

        Input: dict1 is a dict form containing features.

               dict2 is a dict form containing features.

        Output: res is a dict form containing interaction

        features.
    """
    res = {}
    for i in dict1:
        for j in dict2:
            res[i + "*" + j] = round(dict1[i] * dict2[j], 3)
    return res


def CalculateInteraction3(dict1={}, dict2={}):
    """
    Calculate the two interaction features by

    F=[Fa(i)+Fb(i)),Fa(i)*Fb(i)] (2n)

    It's used in same type of descriptors.

    Usage:

        res=CalculateInteraction(dict1,dict2)

        Input: dict1 is a dict form containing features.

               dict2 is a dict form containing features.

        Output: res is a dict form containing interaction

        features.
    """
    res = {}
    for i in dict1:
        res[i + "+" + i] = round(dict1[i] + dict2[i], 3)
        res[i + "*" + i] = round(dict1[i] * dict2[i], 3)
    return res


if __name__ == "__main__":
    import os

    from PyDNA import PyDNAac

    DNA_des = PyDNAac.GetTCC(
        "GACTGAACTGCACTTTGGTTTCATATTATTTGCTC",
        phyche_index=["Dnase I", "Nucleosome", "MW-kg"],
    )

    print(DNA_des)

    from PyProtein import CTD

    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    protein_des = CTD.CalculateCTD(protein)

    from PyMolecule import moe
    from rdkit import Chem

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    m = Chem.MolFromSmiles(smis[3])
    mol_des = moe.GetMOE(m)

    mol_mol_interaction1 = CalculateInteraction1(mol_des, mol_des)
    print(mol_mol_interaction1)

    mol_mol_interaction2 = CalculateInteraction2(mol_des, mol_des)
    print(mol_mol_interaction2)

    mol_mol_interaction3 = CalculateInteraction3(mol_des, mol_des)
    print(mol_mol_interaction3)

    pro_mol_interaction1 = CalculateInteraction1(mol_des, protein_des)
    print(pro_mol_interaction1)

    pro_mol_interaction2 = CalculateInteraction2(mol_des, protein_des)
    print(pro_mol_interaction2)

    DNA_mol_interaction1 = CalculateInteraction1(DNA_des, mol_des)
    print(DNA_mol_interaction1)

    DNA_mol_interaction2 = CalculateInteraction2(DNA_des, mol_des)
    print(DNA_mol_interaction2)
