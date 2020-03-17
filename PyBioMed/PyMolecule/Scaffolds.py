# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
This module provides standard Murcko-type decomposition of molecules into scaffolds.

If you have any question please contact me via email.

Bemis, G. W.; Murcko, M. A. “The Properties of Known Drugs. 1. Molecular Frameworks.”
J. Med. Chem. 39:2887-93 (1996).

2016.11.15

@author: Zhijiang Yao and Dongsheng Cao

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

# Third party modules
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold


def GetScaffold(mol, generic_framework=False):
    """
    #################################################################
    Calculate Scaffold

    Usage:

        result = GetScaffold(mol)

        Input: mol is a molecule object.

        generic_framework is boolean value. If the generic_framework is True, the

        result returns a generic framework.

        Output: result is a string form of the molecule's scaffold.
    #################################################################
    """
    core = MurckoScaffold.GetScaffoldForMol(mol)
    if generic_framework == True:
        fw = MurckoScaffold.MakeScaffoldGeneric(core)
        mol_generic_framework = Chem.MolToSmiles(fw)
        return mol_generic_framework
    else:
        mol_scafflod = Chem.MolToSmiles(core)
        return mol_scafflod


def GetUniqueScaffold(mols, generic_framework=False):
    """
    #################################################################
    Calculate molecules' unique scaffolds

    Usage:

        result = GetUniqueScaffold(mols)

        Input: mols is molecules object.

        generic_framework is boolean value. If the generic_framework is True, the

        result returns a generic framework dictionary.

        Output: result is a tuple form. The first is the list of

        unique scaffolds. The second is a dict form whose keys are the

        scaffolds and the values are the scaffolds' count.
    #################################################################
    """
    scaffolds_dict = {}
    if generic_framework == True:
        for mol in mols:
            scaffold = GetScaffold(mol, generic_framework=True)
            if scaffold not in scaffolds_dict:
                scaffolds_dict[scaffold] = 1
            else:
                scaffolds_dict[scaffold] += 1
    else:
        for mol in mols:
            scaffold = GetScaffold(mol)
            if scaffold not in scaffolds_dict:
                scaffolds_dict[scaffold] = 1
            else:
                scaffolds_dict[scaffold] += 1
    return scaffolds_dict.keys(), scaffolds_dict


if __name__ == "__main__":
    m1 = Chem.MolFromSmiles(
        "O=C1N=C(Nc2ncc(nc12)CNc1ccc(cc1)C(=O)N[C@@H](CCC(O)=O)C(O)=O)N"
    )
    print(GetScaffold(m1))
    print(GetScaffold(m1, generic_framework=True))
    mols = [m1] * 3
    print(GetUniqueScaffold(mols))
