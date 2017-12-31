# -*- coding: utf-8 -*-
"""
###############################################################################

The module is used for computing the composition of amino acids, dipetide and 

3-mers (tri-peptide) for a given protein sequence. You can get 8420 descriptors 

for a given protein sequence. You can freely use and distribute it. If you hava 

any problem, you could contact with us timely!

References:

[1]: Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein

fold class predictions. Nucleic Acids Res, 22, 3616-3619.

[2]: Hua, S. and Sun, Z. (2001) Support vector machine approach for protein

subcellular localization prediction. Bioinformatics, 17, 721-728.


[3]:Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold class

prediction: new methods of statistical classification. Proc Int Conf Intell Syst Mol

Biol, 106-112.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

###############################################################################
"""

import re

AALetter = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


#############################################################################################
def CalculateAAComposition(ProteinSequence):
    """
    ########################################################################
    Calculate the composition of Amino acids

    for a given protein sequence.

    Usage:

    result=CalculateAAComposition(protein)

    Input: protein is a pure protein sequence.

    Output: result is a dict form containing the composition of

    20 amino acids.
    ########################################################################
    """
    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        Result[i] = round(float(ProteinSequence.count(i)) / LengthSequence * 100, 3)
    return Result


#############################################################################################
def CalculateDipeptideComposition(ProteinSequence):
    """
    ########################################################################
    Calculate the composition of dipeptidefor a given protein sequence.

    Usage:

    result=CalculateDipeptideComposition(protein)

    Input: protein is a pure protein sequence.

    Output: result is a dict form containing the composition of

    400 dipeptides.
    ########################################################################
    """

    LengthSequence = len(ProteinSequence)
    Result = {}
    for i in AALetter:
        for j in AALetter:
            Dipeptide = i + j
            Result[Dipeptide] = round(float(ProteinSequence.count(Dipeptide)) / (LengthSequence - 1) * 100, 2)
    return Result


#############################################################################################

def Getkmers():
    """
    ########################################################################
    Get the amino acid list of 3-mers.

    Usage:

    result=Getkmers()

    Output: result is a list form containing 8000 tri-peptides.

    ########################################################################
    """
    kmers = list()
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                kmers.append(i + j + k)
    return kmers


#############################################################################################
def GetSpectrumDict(proteinsequence):
    """
    ########################################################################
    Calcualte the spectrum descriptors of 3-mers for a given protein.

    Usage:

    result=GetSpectrumDict(protein)

    Input: protein is a pure protein sequence.

    Output: result is a dict form containing the composition values of 8000

    3-mers.
    ########################################################################
    """
    result = {}
    kmers = Getkmers()
    for i in kmers:
        result[i] = len(re.findall(i, proteinsequence))
    return result


#############################################################################################
def CalculateAADipeptideComposition(ProteinSequence):
    """
    ########################################################################
    Calculate the composition of AADs, dipeptide and 3-mers for a

    given protein sequence.

    Usage:

    result=CalculateAADipeptideComposition(protein)

    Input: protein is a pure protein sequence.

    Output: result is a dict form containing all composition values of

    AADs, dipeptide and 3-mers (8420).
    ########################################################################
    """

    result = {}
    result.update(CalculateAAComposition(ProteinSequence))
    result.update(CalculateDipeptideComposition(ProteinSequence))
    result.update(GetSpectrumDict(ProteinSequence))

    return result


#############################################################################################
if __name__ == "__main__":
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"

    AAC = CalculateAAComposition(protein)
    print AAC
    DIP = CalculateDipeptideComposition(protein)
    print DIP
    spectrum = GetSpectrumDict(protein)
    print spectrum
    res = CalculateAADipeptideComposition(protein)
    print len(res)
