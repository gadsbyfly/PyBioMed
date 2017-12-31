# -*- coding: utf-8 -*-
"""
The script is used for testing.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com 
"""


import os

from PyBioMed.PyProtein.AAComposition import CalculateAAComposition, CalculateDipeptideComposition, GetSpectrumDict

from PyBioMed.PyProtein.Autocorrelation import CalculateNormalizedMoreauBrotoAutoTotal,CalculateMoranAutoTotal,CalculateGearyAutoTotal

from PyBioMed.PyProtein.Autocorrelation import CalculateEachGearyAuto,CalculateEachMoranAuto,CalculateEachNormalizedMoreauBrotoAuto

from PyBioMed.PyProtein.CTD import CalculateCTD

from PyBioMed.PyProtein.QuasiSequenceOrder import GetSequenceOrderCouplingNumberTotal, GetQuasiSequenceOrder

from PyBioMed.PyProtein.QuasiSequenceOrder import GetSequenceOrderCouplingNumberp,GetQuasiSequenceOrderp

from PyBioMed.PyProtein.PseudoAAC import _GetPseudoAAC, GetAPseudoAAC,GetPseudoAAC

from PyBioMed.PyProtein.GetSubSeq import GetSubSequence

from PyBioMed.PyProtein.AAIndex import GetAAIndex1, GetAAIndex23

from PyBioMed.PyProtein.ConjointTriad import CalculateConjointTriad

from PyBioMed.PyProtein import AAComposition,Autocorrelation,CTD,QuasiSequenceOrder,PseudoAAC,GetProteinFromUniprot,GetSubSeq



modulelists=['AAComposition','Autocorrelation','CTD','QuasiSequenceOrder','PseudoAAC','GetProteinFromUniprot','GetSubSeq']


def test_pyprotein():

    AAC  = eval(modulelists[0])
    AC   = eval(modulelists[1])
    CTD  = eval(modulelists[2])
    QSO  = eval(modulelists[3])
    PAAC = eval(modulelists[4])
    GPFU = eval(modulelists[5])
    GSS  = eval(modulelists[6])
    
   
    
    
    print '...............................................................'
    
    print "testing the GetSubSeq module"
    
    ProteinSequence = 'ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS'
    
    temp=GSS.GetSubSequence(ProteinSequence,ToAA='D', window=5)
    
    print temp
    
    print '...............................................................'
    
    print "testing the AAComposition module"
    
    temp=AAC.CalculateAAComposition(ProteinSequence)
    
    print temp
    
    temp=AAC.CalculateDipeptideComposition(ProteinSequence)
    
    temp=AAC.GetSpectrumDict(ProteinSequence)
    
    temp=AAC.CalculateAADipeptideComposition(ProteinSequence)
    
    print '...............................................................'
    
    print "testing the Autocorrelation module"
    
    temp=AC.CalculateNormalizedMoreauBrotoAuto(ProteinSequence,[AC._ResidueASA],['ResidueASA'])
    
    print temp
    
    temp=AC.CalculateMoranAuto(ProteinSequence,[AC._ResidueASA],['ResidueASA'])
    
    print temp
    
    temp=AC.CalculateGearyAuto(ProteinSequence,[AC._ResidueASA],['ResidueASA'])
    
    print temp
    
    temp=AC.CalculateAutoTotal(ProteinSequence)
    
    print '...............................................................'
    
    print "testing the CTD module"
    
    temp = CTD.CalculateC(ProteinSequence)
    
    print temp
    
    temp = CTD.CalculateT(ProteinSequence)
    
    print temp
    
    temp = CTD.CalculateD(ProteinSequence)
    
    print temp
    
    temp = CTD.CalculateCTD(ProteinSequence)
    
    print temp
    
    print '...............................................................'
    
    print "testing the QuasiSequenceOrder module"
    
    temp=QSO.GetSequenceOrderCouplingNumberTotal(ProteinSequence,maxlag=30)
    
    print temp
    
    temp=QSO.GetQuasiSequenceOrder(ProteinSequence,maxlag=30,weight=0.1)
    
    print temp
    
    print '...............................................................'
    
    print "testing the PseudoAAC module"
    
    temp= PAAC.GetAPseudoAAC(ProteinSequence,lamda=10,weight=0.5)
    
    print temp
    
    
    temp= PAAC._GetPseudoAAC(ProteinSequence,lamda=10,weight=0.05)
    
    print temp
    
    print '...............................................................'
    
    print "Tested successfully!"


if __name__ == '__main__':
    test_pyprotein()






























































































































