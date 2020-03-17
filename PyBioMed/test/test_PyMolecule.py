# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
The script is used for testing.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com
"""
# Core Library modules
import os


def test_pymolecule():

    from rdkit import Chem
    import pybel
    from PyBioMed.PyMolecule.AtomProperty import (
        AtomProperty,
        GetAbsoluteAtomicProperty,
        GetRelativeAtomicProperty,
    )

    from PyBioMed.PyMolecule.estate import GetEstate
    from PyBioMed.PyMolecule.constitution import GetConstitutional
    from PyBioMed.PyMolecule.fingerprint import CalculateECFP4Fingerprint
    from PyBioMed.PyMolecule.connectivity import GetConnectivity
    from PyBioMed.PyMolecule.charge import GetCharge
    from PyBioMed.PyMolecule.bcut import GetBurden
    from PyBioMed.PyMolecule.basak import Getbasak
    from PyBioMed.PyMolecule.geary import GetGearyAuto
    from PyBioMed.PyMolecule.kappa import GetKappa
    from PyBioMed.PyMolecule.moe import GetMOE
    from PyBioMed.PyMolecule.molproperty import GetMolecularProperty
    from PyBioMed.PyMolecule.moran import GetMoranAuto
    from PyBioMed.PyMolecule.moreaubroto import GetMoreauBrotoAuto
    from PyBioMed.PyMolecule.topology import GetTopology
    from PyBioMed.PyMolecule.ghosecrippen import GhoseCrippenFingerprint
    from PyBioMed.PyMolecule.cats2d import CATS2D

    print("...............................................................")
    print("testing the AtomProperty module")

    for i, j in AtomProperty.items():
        print(j)
    print(GetAbsoluteAtomicProperty(element="S", propertyname="En"))
    print(GetRelativeAtomicProperty(element="S", propertyname="En"))
    # ==============================================================================
    # Getbasak
    # ==============================================================================
    print("...............................................................")
    print("testing the Getbasak module")

    smi5 = ["CCCCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCCCCN", "c1ccccc1N"]
    for index, smi in enumerate(smi5):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", Getbasak(m))
        print(len(Getbasak(m)))
    # ==============================================================================
    # GetBurden, bcut model
    # ==============================================================================
    print("...............................................................")
    print("testing the bcut module")

    smi5 = ["CCOCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCCCCN", "c1ccccc1N", "C"]
    for index, smi in enumerate(smi5):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi, "\n")
        print(GetBurden(m))
        print(len(GetBurden(m)))
    # ==============================================================================
    # GetCharge
    # ==============================================================================
    print("...............................................................")
    print("testing the GetCharge module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    smi5 = ["CCCCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCCCCN", "c1ccccc1N"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetCharge(m))
        print(len(GetCharge(m)))
    # ==============================================================================
    # GetConnectivity
    # ==============================================================================
    print("...............................................................")
    print("testing the GetConnectivity module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    smi5 = ["CCCCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCCCCN", "c1ccccc1N"]

    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetConnectivity(m))
        print("\t", len(GetConnectivity(m)))

    # ==============================================================================
    # GetConstitutional
    # ==============================================================================
    print("...............................................................")
    print("testing the GetConstitutional module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    smi5 = ["CCCCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCCCCN", "c1ccccc1N"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetConstitutional(m))
        print(len(GetConstitutional(m)))
    # ==============================================================================
    # GetEstate
    # ==============================================================================
    print("...............................................................")
    print("testing the GetEstate module")

    smi5 = ["COCCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCOCCN", "c1ccccc1N"]
    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        #        print '\t',CalculateEstateFingerprint(m)
        #        print '\t',CalculateEstateValue(m)
        #        print '\t',CalculateMaxAtomTypeEState(m)
        #        print '\t', CalculateMinAtomTypeEState(m)

        print(GetEstate(m))
    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("fingerprint......")

    ms = [
        Chem.MolFromSmiles("CCOC=N"),
        Chem.MolFromSmiles("NC1=NC(=CC=N1)N1C=CC2=C1C=C(O)C=C2"),
    ]
    m2 = [pybel.readstring("smi", "CCOC=N"), pybel.readstring("smi", "CCO")]
    res1 = CalculateECFP4Fingerprint(ms[0])
    print(res1)
    res2 = CalculateECFP4Fingerprint(ms[1])
    print(res2)
    #    print CalculateSimilarityRdkit(res1[2],res2[2])
    #    print CalculateSimilarityPybel(res1,res2)

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetGearyAuto module")

    smi5 = ["COCCCC", "CCC(C)CC", "CC(C)CCC", "CC(C)C(C)C", "CCOCCN", "c1ccccc1N"]
    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        ##        print '\t',CalculateEstateFingerprint(m)
        ##        print '\t',CalculateEstateValue(m)
        ##        print '\t',CalculateMaxAtomTypeEState(m)
        ##        print '\t', CalculateMinAtomTypeEState(m)

        print(GetGearyAuto(m))

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetKappa module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetKappa(m))
        print("\t", len(GetKappa(m)))

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetMOE module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetMOE(m))
        print("\t", len(GetMOE(m)))
    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetMolecularProperty module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetMolecularProperty(m))

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetMoranAuto module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        ##        print '\t',CalculateEstateFingerprint(m)
        ##        print '\t',CalculateEstateValue(m)
        ##        print '\t',CalculateMaxAtomTypeEState(m)
        ##        print '\t', CalculateMinAtomTypeEState(m)

        print(GetMoranAuto(m))

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetMoreauBrotoAuto module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-].[Na+]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        ##        print '\t',CalculateEstateFingerprint(m)
        ##        print '\t',CalculateEstateValue(m)
        ##        print '\t',CalculateMaxAtomTypeEState(m)
        ##        print '\t', CalculateMinAtomTypeEState(m)

        print(len(GetMoreauBrotoAuto(m)))
    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the GetTopology module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetTopology(m))
        print("\t", len(GetTopology(m)))

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the CATS2D module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-]"]

    for i in smis:
        mol = Chem.MolFromSmiles(i)
        cats = CATS2D(mol, PathLength=10, scale=3)
        print("\t", cats)
        print("\t", len(cats))

    # ==============================================================================
    #
    # ==============================================================================
    print("...............................................................")
    print("testing the ghoseFP module")

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-]"]

    for i in smis:
        mol = Chem.MolFromSmiles(i)
        ghoseFP = GhoseCrippenFingerprint(mol)
        print("\t", ghoseFP)
        print("\t", len(ghoseFP))

    print("testing the ghose_count module")

    for i in smis:
        mol = Chem.MolFromSmiles(i)
        ghose_count = GhoseCrippenFingerprint(mol, count=True)
        print("\t", ghose_count)
        print("\t", len(ghose_count))

    # ==============================================================================
    #
    # ==============================================================================


if __name__ == "__main__":
    test_pymolecule()
