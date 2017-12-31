# -*- coding: utf-8 -*-
"""
The script is used for testing.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com 
"""
import os

import string

from PyBioMed.PyGetMol.GetDNA import GetDNAFromUniGene

from PyBioMed.PyGetMol.GetProtein import GetSeqFromPDB
from PyBioMed.PyGetMol.GetProtein import GetPDB

from PyBioMed.PyGetMol.Getmol import ReadMolFromSDF
from PyBioMed.PyGetMol.Getmol import ReadMolFromMOL
from PyBioMed.PyGetMol.Getmol import ReadMolFromSmile
from PyBioMed.PyGetMol.Getmol import GetMolFromCAS
from PyBioMed.PyGetMol.Getmol import GetMolFromNCBI
from PyBioMed.PyGetMol.Getmol import GetMolFromDrugbank
from PyBioMed.PyGetMol.Getmol import GetMolFromKegg
from PyBioMed.PyGetMol.GetProtein import timelimited



def test_pygetmol():
    print '-'*10+'START'+'-'*10
    #==============================================================================
    # GetDNA
    #==============================================================================
    @timelimited(30)
    def test_GetDNAFromUniGene():
        seqid = 'AA954964'
        seqid2 = 'CB216422'

        try:
            sequence1 = GetDNAFromUniGene(seqid)
            sequence2 = GetDNAFromUniGene(seqid2)

            print sequence1
            print sequence2
        except:
            print "Can't visit the internet"
    test_GetDNAFromUniGene()
    print '-'*25
    #==============================================================================
    # Get protein
    #==============================================================================
    @timelimited(30)
    def test_GetPDB():
        try:
            GetPDB(['1atp','1efz','1f88'])

            seq = GetSeqFromPDB('1atp.pdb')
            print seq

            seq2 = GetSeqFromPDB('1efz.pdb')
            print seq2

            seq3 = GetSeqFromPDB('1f88.pdb')
            print seq3
        except:
            print "Can't visit the internet"
    test_GetPDB()
    print '-'*25
    #==============================================================================
    # Get molecule
    #==============================================================================
    @timelimited(30)
    def test_GetSmallMol():
        try:
            temp=GetMolFromCAS(casid="50-12-4")
            print temp
            temp=GetMolFromNCBI(cid="2244")
            print temp
            temp=GetMolFromDrugbank(dbid="DB00133")
            print temp
            temp=GetMolFromKegg(kid="D02176")
            print temp
        except:
            print "Can't visit the internet"
    test_GetSmallMol()

    print '-'*10+'END'+'-'*10

if __name__ == '__main__':
    test_pygetmol()









