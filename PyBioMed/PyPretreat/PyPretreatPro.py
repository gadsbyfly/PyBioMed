# -*- coding: utf-8 -*-
"""
#####################################################################################
This module is used for checking whether the input protein sequence is valid amino acid

sequence. You can freely use and distribute it. If you hava any problem, you could 

contact with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

#####################################################################################

"""

AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


def ProteinCheck(ProteinSequence):
    """
    ###################################################################################
    Check whether the protein sequence is a valid amino acid sequence or not

    Usage:

    result=ProteinCheck(protein)

    Input: protein is a pure protein sequence.

    Output: if the check is no problem, result will return the length of protein.

    if the check has problems, result will return 0.
    ###################################################################################
    """

    NumPro = len(ProteinSequence)
    for i in ProteinSequence:
        if i not in AALetter:
            flag = 0
            break
        else:
            flag = NumPro

    return flag


if __name__ == "__main__":
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDASU"
    protein = "ADGC"
    print ProteinCheck(protein)
