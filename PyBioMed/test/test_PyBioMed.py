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

# First party modules
from PyBioMed.test.test_PyDNA import test_pydna
from PyBioMed.test.test_PyGetMol import test_pygetmol
from PyBioMed.test.test_PyInteration import test_pyinteration
from PyBioMed.test.test_PyMolecule import test_pymolecule
from PyBioMed.test.test_PyPretreat import test_pypretreat
from PyBioMed.test.test_PyProtein import test_pyprotein


def test_pybiomed():

    test_pydna()

    test_pyinteration()

    test_pypretreat()

    test_pygetmol()

    test_pymolecule()

    test_pyprotein()


if __name__ == "__main__":
    test_pybiomed()
