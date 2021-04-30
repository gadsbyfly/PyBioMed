# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
###############################################################################
This module is used for calculating the conjoint triad features only from the

protein sequence information. You can get 7*7*7=343 features.You can freely

use and distribute it. If you hava any problem, you could contact with us timely!

Reference:

Juwen Shen, Jian Zhang, Xiaomin Luo, Weiliang Zhu, Kunqian Yu, Kaixian Chen,

Yixue Li, Huanliang Jiang. Predicting proten-protein interactions based only

on sequences inforamtion. PNAS. 2007 (104) 4337-4341.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

###############################################################################
"""

# Core Library modules
import string

###############################################################################
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

# a Dipole scale (Debye): -, Dipole<1.0; +, 1.0<Dipole<2.0; ++, 2.0<Dipole<3.0; +++, Dipole>3.0; +'+'+', Dipole>3.0 with opposite orientation.
# b Volume scale (Å3): -, Volume<50; +, Volume> 50.
# c Cys is separated from class 3 because of its ability to form disulfide bonds.

_repmat = {
    1: ["A", "G", "V"],
    2: ["I", "L", "F", "P"],
    3: ["Y", "M", "T", "S"],
    4: ["H", "N", "Q", "W"],
    5: ["R", "K"],
    6: ["D", "E"],
    7: ["C"],
}


###############################################################################


def _Str2Num(proteinsequence):
    """
    translate the amino acid letter into the corresponding class based on the

    given form.

    """
    repmat = {}
    for i in _repmat:
        for j in _repmat[i]:
            repmat[j] = i

    res = proteinsequence
    for i in repmat:
        res = res.replace(i, str(repmat[i]))
    return res


###############################################################################
def CalculateConjointTriad(proteinsequence):
    """
    Calculate the conjoint triad features from protein sequence.

    Useage:

    res = CalculateConjointTriad(protein)

    Input: protein is a pure protein sequence.

    Output is a dict form containing all 343 conjoint triad features.
    """
    res = {}
    proteinnum = _Str2Num(proteinsequence)
    for i in range(1, 8):
        for j in range(1, 8):
            for k in range(1, 8):
                temp = str(i) + str(j) + str(k)
                res[temp] = proteinnum.count(temp)
    return res


###############################################################################

if __name__ == "__main__":
    protein = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    print(CalculateConjointTriad(protein))
