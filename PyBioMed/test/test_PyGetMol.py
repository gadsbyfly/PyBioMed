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
import string

# First party modules
from PyBioMed.PyGetMol.GetDNA import GetDNAFromUniGene
from PyBioMed.PyGetMol.Getmol import (
    GetMolFromCAS,
    GetMolFromDrugbank,
    GetMolFromKegg,
    GetMolFromNCBI,
    ReadMolFromMOL,
    ReadMolFromSDF,
    ReadMolFromSmile,
)
from PyBioMed.PyGetMol.GetProtein import GetPDB, GetSeqFromPDB, timelimited


def test_pygetmol():
    print("-" * 10 + "START" + "-" * 10)
    # ==============================================================================
    # GetDNA
    # ==============================================================================
    @timelimited(30)
    def test_GetDNAFromUniGene():
        seqid = "AA954964"
        seqid2 = "CB216422"

        try:
            sequence1 = GetDNAFromUniGene(seqid)
            sequence2 = GetDNAFromUniGene(seqid2)

            print(sequence1)
            print(sequence2)
        except:
            print("Can't visit the internet")

    test_GetDNAFromUniGene()
    print("-" * 25)
    # ==============================================================================
    # Get protein
    # ==============================================================================
    @timelimited(30)
    def test_GetPDB():
        try:
            GetPDB(["1atp", "1efz", "1f88"])

            seq = GetSeqFromPDB("1atp.pdb")
            print(seq)

            seq2 = GetSeqFromPDB("1efz.pdb")
            print(seq2)

            seq3 = GetSeqFromPDB("1f88.pdb")
            print(seq3)
        except:
            print("Can't visit the internet")

    test_GetPDB()
    print("-" * 25)
    # ==============================================================================
    # Get molecule
    # ==============================================================================
    @timelimited(30)
    def test_GetSmallMol():
        try:
            temp = GetMolFromCAS(casid="50-12-4")
            print(temp)
            temp = GetMolFromNCBI(cid="2244")
            print(temp)
            temp = GetMolFromDrugbank(dbid="DB00133")
            print(temp)
            temp = GetMolFromKegg(kid="D02176")
            print(temp)
        except:
            print("Can't visit the internet")

    test_GetSmallMol()

    print("-" * 10 + "END" + "-" * 10)


if __name__ == "__main__":
    test_pygetmol()
