# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
A class used for computing different types of protein descriptors!

You can freely use and distribute it. If you have any problem,

you could contact with us timely.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

"""
# Third party modules
from AAComposition import (
    CalculateAAComposition,
    CalculateDipeptideComposition,
    GetSpectrumDict,
)
from AAIndex import GetAAIndex1, GetAAIndex23
from Autocorrelation import (
    CalculateEachGearyAuto,
    CalculateEachMoranAuto,
    CalculateEachNormalizedMoreauBrotoAuto,
    CalculateGearyAutoTotal,
    CalculateMoranAutoTotal,
    CalculateNormalizedMoreauBrotoAutoTotal,
)
from ConjointTriad import CalculateConjointTriad
from CTD import CalculateCTD
from GetSubSeq import GetSubSequence
from PseudoAAC import GetAPseudoAAC, GetPseudoAAC, _GetPseudoAAC
from QuasiSequenceOrder import (
    GetQuasiSequenceOrder,
    GetQuasiSequenceOrderp,
    GetSequenceOrderCouplingNumberp,
    GetSequenceOrderCouplingNumberTotal,
)


class PyProtein:
    """
    This GetProDes class aims at collecting all descriptor calcualtion modules into a simple class.

    """

    AALetter = [
        "A",
        "R",
        "N",
        "D",
        "C",
        "E",
        "Q",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
    ]

    Version = 1.0

    def __init__(self, ProteinSequence=""):
        """
        input a protein sequence
        """
        if len(ProteinSequence) == 0:
            print(
                "You must input a protein sequence when constructing a object. It is a string!"
            )
        else:
            self.ProteinSequence = ProteinSequence

    def GetAAComp(self):
        """
        amino acid compositon descriptors (20)

        Usage:

        result = GetAAComp()
        """
        res = CalculateAAComposition(self.ProteinSequence)
        return res

    def GetDPComp(self):
        """
        dipeptide composition descriptors (400)

        Usage:

        result = GetDPComp()
        """
        res = CalculateDipeptideComposition(self.ProteinSequence)
        return res

    def GetTPComp(self):
        """
        tri-peptide composition descriptors (8000)

        Usage:

        result = GetTPComp()
        """
        res = GetSpectrumDict(self.ProteinSequence)
        return res

    def GetMoreauBrotoAuto(self):
        """
        Normalized Moreau-Broto autocorrelation descriptors (240)

        Usage:

        result = GetMoreauBrotoAuto()
        """
        res = CalculateNormalizedMoreauBrotoAutoTotal(self.ProteinSequence)
        return res

    def GetMoranAuto(self):
        """
        Moran autocorrelation descriptors (240)

        Usage:

        result = GetMoranAuto()
        """
        res = CalculateMoranAutoTotal(self.ProteinSequence)
        return res

    def GetGearyAuto(self):
        """
        Geary autocorrelation descriptors (240)

        Usage:

        result = GetGearyAuto()
        """
        res = CalculateGearyAutoTotal(self.ProteinSequence)
        return res

    def GetCTD(self):
        """
        Composition Transition Distribution descriptors (147)

        Usage:

        result = GetCTD()
        """
        res = CalculateCTD(self.ProteinSequence)
        return res

    def GetPAAC(self, lamda=10, weight=0.05):
        """
        Type I Pseudo amino acid composition descriptors (default is 30)

        Usage:

        result = GetPAAC(lamda=10,weight=0.05)

        lamda factor reflects the rank of correlation and is a non-Negative integer, such as 15.

        Note that (1)lamda should NOT be larger than the length of input protein sequence;

        (2) lamda must be non-Negative integer, such as 0, 1, 2, ...; (3) when lamda =0, the

        output of PseAA server is the 20-D amino acid composition.

        weight factor is designed for the users to put weight on the additional PseAA components

        with respect to the conventional AA components. The user can select any value within the

        region from 0.05 to 0.7 for the weight factor.
        """
        res = _GetPseudoAAC(self.ProteinSequence, lamda=lamda, weight=weight)
        return res

    def GetPAACp(self, lamda=10, weight=0.05, AAP=[]):
        """
        Type I Pseudo amino acid composition descriptors for the given properties (default is 30)

        Usage:

        result = GetPAACp(lamda=10,weight=0.05,AAP=[])

        lamda factor reflects the rank of correlation and is a non-Negative integer, such as 15.

        Note that (1)lamda should NOT be larger than the length of input protein sequence;

        (2) lamda must be non-Negative integer, such as 0, 1, 2, ...; (3) when lamda =0, the

        output of PseAA server is the 20-D amino acid composition.

        weight factor is designed for the users to put weight on the additional PseAA components

        with respect to the conventional AA components. The user can select any value within the

        region from 0.05 to 0.7 for the weight factor.

        AAP is a list form containing the properties, each of which is a dict form.
        """
        res = GetPseudoAAC(self.ProteinSequence, lamda=lamda, weight=weight, AAP=AAP)
        return res

    def GetAPAAC(self, lamda=10, weight=0.5):
        """
        Amphiphilic (Type II) Pseudo amino acid composition descriptors

        default is 30

        Usage:

        result = GetAPAAC(lamda=10,weight=0.5)

        lamda factor reflects the rank of correlation and is a non-Negative integer, such as 15.

        Note that (1)lamda should NOT be larger than the length of input protein sequence;

        (2) lamda must be non-Negative integer, such as 0, 1, 2, ...; (3) when lamda =0, the

        output of PseAA server is the 20-D amino acid composition.

        weight factor is designed for the users to put weight on the additional PseAA components

        with respect to the conventional AA components. The user can select any value within the

        region from 0.05 to 0.7 for the weight factor.

        """
        res = GetAPseudoAAC(self.ProteinSequence, lamda=lamda, weight=weight)
        return res

    def GetSOCN(self, maxlag=45):
        """
        Sequence order coupling numbers  default is 45

        Usage:

        result = GetSOCN(maxlag=45)

        maxlag is the maximum lag and the length of the protein should be larger

        than maxlag. default is 45.
        """
        res = GetSequenceOrderCouplingNumberTotal(self.ProteinSequence, maxlag=maxlag)
        return res

    def GetSOCNp(self, maxlag=45, distancematrix={}):
        """
        Sequence order coupling numbers  default is 45

        Usage:

        result = GetSOCN(maxlag=45)

        maxlag is the maximum lag and the length of the protein should be larger

        than maxlag. default is 45.

        distancematrix is a dict form containing 400 distance values
        """
        res = GetSequenceOrderCouplingNumberp(
            self.ProteinSequence, maxlag=maxlag, distancematrix=distancematrix
        )
        return res

    def GetQSO(self, maxlag=30, weight=0.1):
        """
        Quasi sequence order descriptors  default is 50

        result = GetQSO(maxlag=30, weight=0.1)

        maxlag is the maximum lag and the length of the protein should be larger

        than maxlag. default is 45.
        """
        res = GetQuasiSequenceOrder(self.ProteinSequence, maxlag=maxlag, weight=weight)
        return res

    def GetQSOp(self, maxlag=30, weight=0.1, distancematrix={}):
        """
        Quasi sequence order descriptors  default is 50

        result = GetQSO(maxlag=30, weight=0.1)

        maxlag is the maximum lag and the length of the protein should be larger

        than maxlag. default is 45.

        distancematrix is a dict form containing 400 distance values
        """
        res = GetQuasiSequenceOrderp(
            self.ProteinSequence,
            maxlag=maxlag,
            weight=weight,
            distancematrix=distancematrix,
        )
        return res

    def GetMoreauBrotoAutop(self, AAP={}, AAPName="p"):
        """
        Normalized Moreau-Broto autocorrelation descriptors for the given property (30)

        Usage:

        result = GetMoreauBrotoAutop(AAP={},AAPName='p')

        AAP is a dict containing physicochemical properities of 20 amino acids
        """
        res = CalculateEachNormalizedMoreauBrotoAuto(
            self.ProteinSequence, AAP=AAP, AAPName=AAPName
        )
        return res

    def GetMoranAutop(self, AAP={}, AAPName="p"):
        """
        Moran autocorrelation descriptors for the given property (30)

        Usage:

        result = GetMoranAutop(AAP={},AAPName='p')

        AAP is a dict containing physicochemical properities of 20 amino acids
        """
        res = CalculateEachMoranAuto(self.ProteinSequence, AAP=AAP, AAPName=AAPName)
        return res

    def GetGearyAutop(self, AAP={}, AAPName="p"):
        """
        Geary autocorrelation descriptors for the given property (30)

        Usage:

        result = GetGearyAutop(AAP={},AAPName='p')

        AAP is a dict containing physicochemical properities of 20 amino acids
        """
        res = CalculateEachGearyAuto(self.ProteinSequence, AAP=AAP, AAPName=AAPName)
        return res

    def GetSubSeq(self, ToAA="S", window=3):
        """
        obtain the sub sequences wit length 2*window+1, whose central point is ToAA

        Usage:

        result = GetSubSeq(ToAA='S',window=3)

        ToAA is the central (query point) amino acid in the sub-sequence.

        window is the span.
        """
        res = GetSubSequence(self.ProteinSequence, ToAA=ToAA, window=window)
        return res

    def GetTriad(self):
        """
        Calculate the conjoint triad features from protein sequence.

        Useage:

        res = GetTriad()

        Output is a dict form containing all 343 conjoint triad features.
        """
        res = CalculateConjointTriad(self.ProteinSequence)
        return res

    def GetALL(self):
        """
        Calcualte all descriptors except tri-peptide descriptors
        """
        res = {}
        res.update(self.GetAAComp())
        res.update(self.GetDPComp())
        res.update(self.GetTPComp())
        res.update(self.GetMoreauBrotoAuto())
        res.update(self.GetMoranAuto())
        res.update(self.GetGearyAuto())
        res.update(self.GetCTD())
        res.update(self.GetPAAC())
        res.update(self.GetAPAAC())
        res.update(self.GetSOCN())
        res.update(self.GetQSO())
        res.update(self.GetTriad())
        return res

    def GetAAindex1(self, name, path="."):
        """
        Get the amino acid property values from aaindex1

        Usage:

        result=GetAAIndex1(name)

        Input: name is the name of amino acid property (e.g., KRIW790103)

        Output: result is a dict form containing the properties of 20 amino acids
        """

        return GetAAIndex1(name, path=path)

    def GetAAindex23(self, name, path="."):
        """
        Get the amino acid property values from aaindex2 and aaindex3

        Usage:

        result=GetAAIndex23(name)

        Input: name is the name of amino acid property (e.g.,TANS760101,GRAR740104)

        Output: result is a dict form containing the properties of 400 amino acid pairs
        """
        return GetAAIndex23(name, path=path)


#####################################################################################################
if __name__ == "__main__":
    import os
    from Autocorrelation import _Steric
    from PseudoAAC import _Hydrophobicity, _hydrophilicity, _residuemass

    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    cds = PyProtein(protein)

    print(cds.GetAAComp())
    print(cds.GetDPComp())
    print(cds.GetTPComp())
    print(cds.GetCTD())
    print(cds.GetPAAC(lamda=5))
    print(cds.GetALL())
    print(cds.GetMoreauBrotoAutop(AAP=_Steric, AAPName="Steric"))
    print(cds.GetMoranAutop(AAP=_Steric, AAPName="Steric"))
    print(cds.GetGearyAutop(AAP=_Steric, AAPName="Steric"))

    print(cds.GetPAACp(lamda=5, weight=0.05, AAP=[_Hydrophobicity, _hydrophilicity]))

    print(cds.GetSubSeq(ToAA="D", window=5))
    print(cds.GetTriad())
    proper = cds.GetAAindex23("GRAR740104")
    print(cds.GetAAindex1("KRIW790103"))

    print(cds.GetQSOp(maxlag=30, weight=0.1, distancematrix=proper))
    print(cds.GetSOCNp(maxlag=30, distancematrix=proper))
