# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
The calculation of Kier and Hall's kappa indices based on its topological

structure. You can get 7 molecular kappa descriptors. You can

freely use and distribute it. If you hava  any problem, you could contact

with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""


# Third party modules
from rdkit import Chem
from rdkit.Chem import pyPeriodicTable as PeriodicTable
from rdkit.Chem import rdchem

periodicTable = rdchem.GetPeriodicTable()


Version = 1.0
################################################################


def CalculateKappa1(mol):
    """
    #################################################################
    Calculation of molecular shape index for one bonded fragment

    ---->kappa1

    Usage:

        result=CalculateKappa1(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    P1 = mol.GetNumBonds(onlyHeavy=1)
    A = mol.GetNumAtoms(onlyHeavy=1)
    denom = P1 + 0.0
    if denom:
        kappa = (A) * (A - 1) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)


def CalculateKappa2(mol):
    """
    #################################################################
    Calculation of molecular shape index for two bonded fragment

    ---->kappa2

    Usage:

        result=CalculateKappa2(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
    A = mol.GetNumAtoms(onlyHeavy=1)

    denom = P2 + 0.0
    if denom:
        kappa = (A - 1) * (A - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)


def CalculateKappa3(mol):
    """
    #################################################################
    Calculation of molecular shape index for three bonded fragment

    ---->kappa3

    Usage:

        result=CalculateKappa3(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
    A = mol.GetNumAtoms(onlyHeavy=1)

    denom = P3 + 0.0
    if denom:
        if A % 2 == 1:
            kappa = (A - 1) * (A - 3) ** 2 / denom ** 2
        else:
            kappa = (A - 3) * (A - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)


def _HallKierAlpha(mol):
    """
    #################################################################
    *Internal Use Only*

    Calculation of the Hall-Kier alpha value for a molecule
    #################################################################
    """
    alphaSum = 0.0
    rC = PeriodicTable.nameTable["C"][5]
    for atom in mol.GetAtoms():
        atNum = atom.GetAtomicNum()
        if not atNum:
            continue
        symb = atom.GetSymbol()
        alphaV = PeriodicTable.hallKierAlphas.get(symb, None)
        if alphaV is not None:
            hyb = atom.GetHybridization() - 2
            if hyb < len(alphaV):
                alpha = alphaV[hyb]
                if alpha is None:
                    alpha = alphaV[-1]
            else:
                alpha = alphaV[-1]
        else:
            rA = PeriodicTable.nameTable[symb][5]
            alpha = rA / rC - 1
        alphaSum += alpha
    return alphaSum


def CalculateKappaAlapha1(mol):
    """
    #################################################################
    Calculation of molecular shape index for one bonded fragment

    with Alapha

    ---->kappam1

    Usage:

        result=CalculateKappaAlapha1(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    P1 = mol.GetNumBonds(onlyHeavy=1)
    A = mol.GetNumAtoms(onlyHeavy=1)
    alpha = _HallKierAlpha(mol)
    denom = P1 + alpha
    if denom:
        kappa = (A + alpha) * (A + alpha - 1) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)


def CalculateKappaAlapha2(mol):
    """
    #################################################################
    Calculation of molecular shape index for two bonded fragment

    with Alapha

    ---->kappam2

    Usage:

        result=CalculateKappaAlapha2(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
    A = mol.GetNumAtoms(onlyHeavy=1)
    alpha = _HallKierAlpha(mol)
    denom = P2 + alpha
    if denom:
        kappa = (A + alpha - 1) * (A + alpha - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)


def CalculateKappaAlapha3(mol):
    """
    #################################################################
    Calculation of molecular shape index for three bonded fragment

    with Alapha

    ---->kappam3

    Usage:

        result=CalculateKappaAlapha3(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
    A = mol.GetNumAtoms(onlyHeavy=1)
    alpha = _HallKierAlpha(mol)
    denom = P3 + alpha
    if denom:
        if A % 2 == 1:
            kappa = (A + alpha - 1) * (A + alpha - 3) ** 2 / denom ** 2
        else:
            kappa = (A + alpha - 3) * (A + alpha - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)


def CalculateFlexibility(mol):
    """
    #################################################################
    Calculation of Kier molecular flexibility index

    ---->phi

    Usage:

        result=CalculateFlexibility(mol)

        Input: mol is a molecule object.

        Output: result is a numeric value.
    #################################################################
    """
    kappa1 = CalculateKappaAlapha1(mol)
    kappa2 = CalculateKappaAlapha2(mol)
    A = mol.GetNumAtoms(onlyHeavy=1)
    phi = kappa1 * kappa2 / (A + 0.0)
    return phi


def GetKappa(mol):
    """
    #################################################################
    Calculation of all kappa values.

    Usage:

        result=GetKappa(mol)

        Input: mol is a molecule object.

        Output: result is a dcit form containing 6 kappa values.
    #################################################################
    """
    res = {}
    res["kappa1"] = CalculateKappa1(mol)
    res["kappa2"] = CalculateKappa2(mol)
    res["kappa3"] = CalculateKappa3(mol)
    res["kappam1"] = CalculateKappaAlapha1(mol)
    res["kappam2"] = CalculateKappaAlapha2(mol)
    res["kappam3"] = CalculateKappaAlapha3(mol)
    res["phi"] = CalculateFlexibility(mol)
    return res


################################################################
if __name__ == "__main__":

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetKappa(m))
        print("\t", len(GetKappa(m)))
