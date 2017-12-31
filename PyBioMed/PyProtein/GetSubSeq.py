# -*- coding: utf-8 -*-
"""
#####################################################################################

The prediction of functional sites (e.g.,methylation) of proteins usually needs to 

split the total protein into a set of segments around specific amino acid. Given a 

specific window size p, we can obtain all segments of length equal to (2*p+1) very 

easily. Note that the output of the method is a list form. You can freely use and 

distribute it. If you have any problem, you could contact with us timely.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

#####################################################################################

"""

import re
import string

AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


#############################################################################################

def GetSubSequence(ProteinSequence, ToAA='S', window=3):
    """
    #######################################################################
    Get all 2*window+1 sub-sequences whose cener is ToAA in a protein.

    Usage:

    result=GetSubSequence(protein,ToAA,window)

    Input:protein is a pure problem sequence.

    ToAA is the central (query point) amino acid in the sub-sequence.

    window is the span.

    result is a list form containing all satisfied sub-sequences.
    #######################################################################
    """

    if ToAA not in AALetter:
        ToAA = ProteinSequence[1]

    Num = len(ProteinSequence)
    seqiter = re.finditer(ToAA, ProteinSequence)
    AAindex = []
    for i in seqiter:
        AAindex.append(i.end())

    result = []
    for i in AAindex:
        if i - window > 0 and Num - i + 1 - window > 0:
            temp = ProteinSequence[i - window - 1:i + window]
            result.append(temp)

    return result


#############################################################################################
if __name__ == "__main__":
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    subseq = GetSubSequence(protein, ToAA='D', window=10)
    print subseq
    print len(subseq)
