# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################

A class used for computing different types of DNA descriptors!

You can freely use and distribute it. If you have any problem,

you could contact with us timely.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

# Core Library modules
import os
import pickle
import sys
from math import pow

# First party modules
from PyBioMed.PyDNA.PyDNAnacutil import MakeKmerList
from PyBioMed.PyDNA.PyDNAutil import Frequency

ALPHABET = "ACGT"


def ExtendPhycheIndex(original_index, extend_index):
    """Extend {phyche:[value, ... ]}"""
    if extend_index is None or len(extend_index) == 0:
        return original_index
    for key in list(original_index.keys()):
        original_index[key].extend(extend_index[key])
    return original_index


global phyche_factor_dic1
global phyche_factor_dic2
full_path = os.path.realpath(__file__)

with open("%s/mmc3.data" % os.path.dirname(full_path), "rb") as f:
    phyche_factor_dic1 = pickle.load(f)

with open("%s/mmc4.data" % os.path.dirname(full_path), "rb") as f:
    phyche_factor_dic2 = pickle.load(f)


def GetPhycheFactorDic(k):
    """Get all {nucleotide: [(phyche, value), ...]} dict."""
    global phyche_factor_dic1
    global phyche_factor_dic2
    if 2 == k:
        phyche_factor_dic = phyche_factor_dic1
    elif 3 == k:
        phyche_factor_dic = phyche_factor_dic2
    else:
        sys.stderr.write("The k can just be 2 or 3.")
        sys.exit(0)

    return phyche_factor_dic


def GetPhycheIndex(k, phyche_list):
    """get phyche_value according phyche_list."""
    phyche_value = {}
    if 0 == len(phyche_list):
        for nucleotide in MakeKmerList(k, ALPHABET):
            phyche_value[nucleotide] = []
        return phyche_value

    nucleotide_phyche_value = GetPhycheFactorDic(k)
    for nucleotide in MakeKmerList(k, ALPHABET):
        if nucleotide not in phyche_value:
            phyche_value[nucleotide] = []
        for e in nucleotide_phyche_value[nucleotide]:
            if e[0] in phyche_list:
                phyche_value[nucleotide].append(e[1])

    return phyche_value


def ParallelCorFunction(nucleotide1, nucleotide2, phyche_index):
    """Get the cFactor.(Type1)"""
    temp_sum = 0.0
    phyche_index_values = list(phyche_index.values())
    len_phyche_index = len(phyche_index_values[0])
    for u in range(len_phyche_index):
        temp_sum += pow(
            float(phyche_index[nucleotide1][u]) - float(phyche_index[nucleotide2][u]), 2
        )

    return temp_sum / len_phyche_index


def SeriesCorFunction(nucleotide1, nucleotide2, big_lamada, phyche_value):
    """Get the series correlation Factor(Type 2)."""
    return float(phyche_value[nucleotide1][big_lamada]) * float(
        phyche_value[nucleotide2][big_lamada]
    )


def GetParallelFactor(k, lamada, sequence, phyche_value):
    """Get the corresponding factor theta list."""
    theta = []
    l = len(sequence)

    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - k - i + 1):
            nucleotide1 = sequence[j : j + k]
            nucleotide2 = sequence[j + i : j + i + k]
            temp_sum += ParallelCorFunction(nucleotide1, nucleotide2, phyche_value)

        theta.append(temp_sum / (l - k - i + 1))

    return theta


def GetSeriesFactor(k, lamada, sequence, phyche_value):
    """Get the corresponding series factor theta list."""
    theta = []
    l_seq = len(sequence)
    temp_values = list(phyche_value.values())
    max_big_lamada = len(temp_values[0])

    for small_lamada in range(1, lamada + 1):
        for big_lamada in range(max_big_lamada):
            temp_sum = 0.0
            for i in range(0, l_seq - k - small_lamada + 1):
                nucleotide1 = sequence[i : i + k]
                nucleotide2 = sequence[i + small_lamada : i + small_lamada + k]
                temp_sum += SeriesCorFunction(
                    nucleotide1, nucleotide2, big_lamada, phyche_value
                )

            theta.append(temp_sum / (l_seq - k - small_lamada + 1))

    return theta


def MakePsekncVector(sequence_list, lamada, w, k, phyche_value, theta_type=1):
    """Generate the pseknc vector."""
    kmer = MakeKmerList(k, ALPHABET)
    vector = []

    for sequence in sequence_list:
        if len(sequence) < k or lamada + k > len(sequence):
            error_info = "Sorry, the sequence length must be larger than " + str(
                lamada + k
            )
            sys.stderr.write(error_info)
            sys.exit(0)

        # Get the nucleotide frequency in the DNA sequence.
        fre_list = [Frequency(sequence, str(key)) for key in kmer]
        fre_sum = float(sum(fre_list))

        # Get the normalized occurrence frequency of nucleotide in the DNA sequence.
        fre_list = [e / fre_sum for e in fre_list]

        # Get the theta_list according the Equation 5.
        if 1 == theta_type:
            theta_list = GetParallelFactor(k, lamada, sequence, phyche_value)
        elif 2 == theta_type:
            theta_list = GetSeriesFactor(k, lamada, sequence, phyche_value)
        theta_sum = sum(theta_list)

        # Generate the vector according the Equation 9.
        denominator = 1 + w * theta_sum

        temp_vec = [round(f / denominator, 3) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 4))

        vector.append(temp_vec)

    return vector


def GetParallelFactorPsednc(lamada, sequence, phyche_value):
    """Get the corresponding factor theta list.
       This def is just for dinucleotide."""
    theta = []
    l = len(sequence)

    for i in range(1, lamada + 1):
        temp_sum = 0.0
        for j in range(0, l - 1 - lamada):
            nucleotide1 = sequence[j] + sequence[j + 1]
            nucleotide2 = sequence[j + i] + sequence[j + i + 1]
            temp_sum += ParallelCorFunction(nucleotide1, nucleotide2, phyche_value)

        theta.append(temp_sum / (l - i - 1))

    return theta


def MakeOldPsekncVector(sequence_list, lamada, w, k, phyche_value, theta_type=1):
    """Generate the pseknc vector."""
    kmer = MakeKmerList(k, ALPHABET)
    vector = []

    for sequence in sequence_list:
        if len(sequence) < k or lamada + k > len(sequence):
            error_info = "Sorry, the sequence length must be larger than " + str(
                lamada + k
            )
            sys.stderr.write(error_info)
            sys.exit(0)

        # Get the nucleotide frequency in the DNA sequence.
        fre_list = [Frequency(sequence, str(key)) for key in kmer]
        fre_sum = float(sum(fre_list))

        # Get the normalized occurrence frequency of nucleotide in the DNA sequence.
        fre_list = [e / fre_sum for e in fre_list]

        # Get the theta_list according the Equation 5.
        if 1 == theta_type:
            theta_list = GetParallelFactorPsednc(lamada, sequence, phyche_value)
        elif 2 == theta_type:
            theta_list = GetSeriesFactor(k, lamada, sequence, phyche_value)
        theta_sum = sum(theta_list)

        # Generate the vector according the Equation 9.
        denominator = 1 + w * theta_sum

        temp_vec = [round(f / denominator, 3) for f in fre_list]
        for theta in theta_list:
            temp_vec.append(round(w * theta / denominator, 4))

        vector.append(temp_vec)

    return vector


if __name__ == "__main__":
    # get_phyche_index(2, ['Base stacking'])
    extra_phyche_index = {
        "AA": [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1],
        "AC": [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
        "AG": [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
        "AT": [1.07, 0.22, 0.62, -1.02, 2.51, 1.17, 1],
        "CA": [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
        "CC": [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
        "CG": [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39, 1],
        "CT": [0.78, 0.36, 0.09, 0.68, -0.24, -0.62, 1],
        "GA": [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
        "GC": [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59, 1],
        "GG": [0.06, 1.08, 0.09, 0.56, -0.82, 0.24, 1],
        "GT": [1.50, 0.50, 0.80, 0.13, 1.29, 1.04, 1],
        "TA": [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39, 1],
        "TC": [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71, 1],
        "TG": [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25, 1],
        "TT": [0.06, 0.5, 0.27, 1.59, 0.11, -0.11, 1],
    }
    phyche_index = ExtendPhycheIndex(
        GetPhycheIndex(k=2, phyche_list=["Base stacking", "DNA denaturation"]),
        extra_phyche_index,
    )
    # print(phyche_index)
    import pandas as pd

    df = pd.DataFrame(GetPhycheFactorDic(3))
