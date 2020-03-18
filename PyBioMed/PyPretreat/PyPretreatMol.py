# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
# from PyPretreatMolutil import Standardizer
# from PyPretreatMolutil import validate_smiles
# from PyPretreatMolutil import standardize_smiles
# from PyPretreatMolutil import Validator
# Core Library modules
import logging

# Third party modules
from rdkit import Chem

# First party modules
from PyBioMed.PyPretreat.PyPretreatMolutil import *

log = logging.getLogger(__name__)

map_dict = {
    "1": "disconnect_metals",
    "2": "normalize",
    "3": "addhs",
    "4": "rmhs",
    "5": "reionize",
    "6": "uncharge",
    "7": "largest_fragment",
    "8": "canonicalize_tautomer",
}

# NORMALIZATIONS = NORMALIZATIONS


class StandardizeMol(object):
    """
    The main class for performing standardization of molecules and deriving parent molecules.

    The primary usage is via the :meth:`~molvs.standardize.Standardizer.standardize` method::

        s = Standardizer()
        mol1 = Chem.MolFromSmiles('C1=CC=CC=C1')
        mol2 = s.standardize(mol1)

    There are separate methods to derive fragment, charge, tautomer, isotope and stereo parent molecules.

    """

    def __init__(
        self,
        normalizations=NORMALIZATIONS,
        acid_base_pairs=ACID_BASE_PAIRS,
        tautomer_transforms=TAUTOMER_TRANSFORMS,
        tautomer_scores=TAUTOMER_SCORES,
        max_restarts=MAX_RESTARTS,
        max_tautomers=MAX_TAUTOMERS,
        prefer_organic=PREFER_ORGANIC,
    ):
        """Initialize a Standardizer with optional custom parameters.

        :param normalizations: A list of Normalizations to apply (default: :data:`~molvs.normalize.NORMALIZATIONS`).
        :param acid_base_pairs: A list of AcidBasePairs for competitive reionization (default:
                                :data:`~molvs.charge.ACID_BASE_PAIRS`).
        :param tautomer_transforms: A list of TautomerTransforms to apply (default:
                                    :data:`~molvs.tautomer.TAUTOMER_TRANSFORMS`).
        :param tautomer_scores: A list of TautomerScores used to determine canonical tautomer (default:
                                :data:`~molvs.tautomer.TAUTOMER_SCORES`).
        :param max_restarts: The maximum number of times to attempt to apply the series of normalizations (default 200).
        :param max_tautomers: The maximum number of tautomers to enumerate (default 1000).
        :param prefer_organic: Whether to prioritize organic fragments when choosing fragment parent (default False).
        """
        log.debug("Initializing Standardizer")
        self.normalizations = normalizations
        self.acid_base_pairs = acid_base_pairs
        self.tautomer_transforms = tautomer_transforms
        self.tautomer_scores = tautomer_scores
        self.max_restarts = max_restarts
        self.max_tautomers = max_tautomers
        self.prefer_organic = prefer_organic

    def __call__(self, mol):
        """Calling a Standardizer instance like a function is the same as calling its
        :meth:`~molvs.standardize.Standardizer.standardize` method."""
        return self.standardize(mol)

    def addhs(self, mol):
        from rdkit.Chem import AddHs

        return AddHs(mol)

    def rmhs(self, mol):
        from rdkit.Chem import RemoveHs

        return RemoveHs(mol)

    @memoized_property
    def disconnect_metals(self):
        """
        :returns: A callable :class:`~molvs.metal.MetalDisconnector` instance.
        """
        return MetalDisconnector()

    @memoized_property
    def normalize(self):
        """
        :returns: A callable :class:`~molvs.normalize.Normalizer` instance.
        """
        return Normalizer(
            normalizations=self.normalizations, max_restarts=self.max_restarts
        )

    @memoized_property
    def reionize(self):
        """
        :returns: A callable :class:`~molvs.charge.Reionizer` instance.
        """
        return Reionizer(acid_base_pairs=self.acid_base_pairs)

    @memoized_property
    def uncharge(self):
        """
        :returns: A callable :class:`~molvs.charge.Uncharger` instance.
        """
        return Uncharger()

    @memoized_property
    def largest_fragment(self):
        """
        :returns: A callable :class:`~molvs.fragment.LargestFragmentChooser` instance.
        """
        return LargestFragmentChooser(prefer_organic=self.prefer_organic)

    @memoized_property
    def canonicalize_tautomer(self):
        """
        :returns: A callable :class:`~molvs.tautomer.TautomerCanonicalizer` instance.
        """
        return TautomerCanonicalizer(
            transforms=self.tautomer_transforms,
            scores=self.tautomer_scores,
            max_tautomers=self.max_tautomers,
        )


def StandardMol(mol):
    """
    The function for performing standardization of molecules and deriving parent molecules.
    The function contains derive fragment, charge, tautomer, isotope and stereo parent molecules.
    The primary usage is::

        mol1 = Chem.MolFromSmiles('C1=CC=CC=C1')
        mol2 = s.standardize(mol1)

    """
    s = Standardizer()
    mol = s.disconnect_metals(mol)
    mol = s.normalize(mol)
    mol = s.uncharge(mol)
    mol = s.largest_fragment(mol)
    mol = s.canonicalize_tautomer(mol)
    mol = s.reionize(mol)
    mol = s.addhs(mol)
    mol = s.rmhs(mol)
    return mol


def StandardSmi(smi):
    """
    The function for performing standardization of molecules and deriving parent molecules.
    The function contains derive fragment, charge, tautomer, isotope and stereo parent molecules.
    The primary usage is::

        smi = StandardSmi('C[n+]1c([N-](C))cccc1')

    """
    mol = Chem.MolFromSmiles(smi)
    mol = StandardMol(mol)
    smi = Chem.MolToSmiles(mol, isomericSmiles=True)
    return smi


def ValidatorMol(mol):
    """
    Return log messages for a given SMILES string using the default validations.

    Note: This is a convenience function for quickly validating a single SMILES string.

    :param string smiles: The SMILES for the molecule.
    :returns: A list of log messages.
    :rtype: list of strings.

    """
    return Validator().validate(mol)


def ValidatorSmi(smi):
    """
    Return log messages for a given SMILES string using the default validations.

    Note: This is a convenience function for quickly validating a single SMILES string.

    :param string smiles: The SMILES for the molecule.
    :returns: A list of log messages.
    :rtype: list of strings.

    """

    return validate_smiles(smi)


if __name__ == "__main__":
    smiles = ["O=C([O-])c1ccccc1", "C[n+]1c([N-](C))cccc1", "[2H]C(Cl)(Cl)Cl"]
    mol = Chem.MolFromSmiles("[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1")
    sm = StandardizeMol()
    mol = sm.addhs(mol)
    mol = sm.disconnect_metals(mol)
    mol = sm.largest_fragment(mol)
    mol = sm.normalize(mol)
    mol = sm.uncharge(mol)
    mol = sm.canonicalize_tautomer(mol)
    mol = sm.reionize(mol)
    mol = sm.rmhs(mol)
    mol = sm.addhs(mol)
    print(Chem.MolToSmiles(mol, isomericSmiles=True))
