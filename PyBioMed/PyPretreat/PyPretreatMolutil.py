# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 10:13:55 2016

@author: yzj
molvs.tautomer
~~~~~~~~~~~~~~

This module contains tools for enumerating tautomers and determining a canonical tautomer.

:copyright: Copyright 2014 by Matt Swain.
:license: MIT, see LICENSE file for more details.
"""
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
import copy
import logging
import functools
from itertools import izip, tee
import sys

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondDir, BondStereo, BondType

log = logging.getLogger(__name__)
#==============================================================================
# from .utils import memoized_property
#==============================================================================
def memoized_property(fget):
    """Decorator to create memoized properties."""
    attr_name = '_{}'.format(fget.__name__)

    @functools.wraps(fget)
    def fget_memoized(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fget(self))
        return getattr(self, attr_name)
    return property(fget_memoized)

def pairwise(iterable):
    """Utility function to iterate in a pairwise fashion."""
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)
#==============================================================================
# from .metal import MetalDisconnector
#==============================================================================
class MetalDisconnector(object):
    """Class for breaking covalent bonds between metals and organic atoms under certain conditions."""

    def __init__(self):
        log.debug('Initializing MetalDisconnector')
        # Initialize SMARTS to identify relevant substructures
        # TODO: Use atomic numbers instead of element symbols in SMARTS to allow for isotopes?
        self._metal_nof = Chem.MolFromSmarts('[Li,Na,K,Rb,Cs,F,Be,Mg,Ca,Sr,Ba,Ra,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Al,Ga,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg,Tl,Pb,Bi]~[N,O,F]'.encode('utf8'))
        self._metal_non = Chem.MolFromSmarts('[Al,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,Hf,Ta,W,Re,Os,Ir,Pt,Au]~[B,C,Si,P,As,Sb,S,Se,Te,Cl,Br,I,At]'.encode('utf8'))
        self._free_metal = Chem.MolFromSmarts('[Li,Na,K,Mg,CaX0+0]'.encode('utf8'))
        self._carboxylic = Chem.MolFromSmarts('[CX3](=O)[OX2H1]'.encode('utf8'))

    def __call__(self, mol):
        """Calling a MetalDisconnector instance like a function is the same as calling its disconnect(mol) method."""
        return self.disconnect(mol)

    def disconnect(self, mol):
        """Break covalent bonds between metals and organic atoms under certain conditions.

        The algorithm works as follows:

        - Disconnect N, O, F from any metal.
        - Disconnect other non-metals from transition metals + Al (but not Hg, Ga, Ge, In, Sn, As, Tl, Pb, Bi, Po).
        - For every bond broken, adjust the charges of the begin and end atoms accordingly.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The molecule with metals disconnected.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running MetalDisconnector')
        # TODO: In future, maybe insert zero-order complex/ionic/dative bonds instead of disconnecting?
        # Remove bonds that match SMARTS
        for smarts in [self._metal_nof, self._metal_non]:
            pairs = mol.GetSubstructMatches(smarts)
            edmol = Chem.EditableMol(mol)
            orders = []
            for i, j in pairs:
                # TODO: Could get the valence contributions of the bond instead of GetBondTypeAsDouble?
                orders.append(int(mol.GetBondBetweenAtoms(i, j).GetBondTypeAsDouble()))
                edmol.RemoveBond(i, j)
            # Adjust neighbouring charges accordingly
            mol = edmol.GetMol()
            for n, (i, j) in enumerate(pairs):
                chg = orders[n]
                atom1 = mol.GetAtomWithIdx(i)
                atom1.SetFormalCharge(atom1.GetFormalCharge() + chg)
                atom2 = mol.GetAtomWithIdx(j)
                atom2.SetFormalCharge(atom2.GetFormalCharge() - chg)
                log.info('Removed covalent bond between %s and %s', atom1.GetSymbol(), atom2.GetSymbol())
        # Ionize a free neutral metal with a protonated carboxylic acid
        # TODO: Extend this to other acids?
        # TODO: Move to charge module?
        fms = mol.GetSubstructMatches(self._free_metal)
        carbs = mol.GetSubstructMatches(self._carboxylic)
        if len(fms) == len(carbs) > 0:
            for fm in fms:
                atom = mol.GetAtomWithIdx(fm[0])
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
            for carb in carbs:
                atom = mol.GetAtomWithIdx(carb[2])
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
        log.debug(Chem.MolToSmiles(mol))
        Chem.SanitizeMol(mol)
        return mol
#==============================================================================
#from .fragment import PREFER_ORGANIC, LargestFragmentChooser, FragmentRemover
#==============================================================================
class FragmentPattern(object):
    """A fragment defined by a SMARTS pattern."""

    def __init__(self, name, smarts):
        """Initialize a FragmentPattern with a name and a SMARTS pattern.

        :param name: A name for this FragmentPattern.
        :param smarts: A SMARTS pattern.
        """
        self.name = name
        self.smarts_str = smarts

    @memoized_property
    def smarts(self):
        return Chem.MolFromSmarts(self.smarts_str.encode('utf8'))

    def __repr__(self):
        return 'FragmentPattern({!r}, {!r})'.format(self.name, self.smarts_str)

    def __str__(self):
        return self.name


#: The default list of :class:`FragmentPatterns <molvs.fragment.FragmentPattern>` to be used by
#: :class:`~molvs.fragment.FragmentRemover`.
REMOVE_FRAGMENTS = (
    FragmentPattern('fluorine', '[F]'),
    FragmentPattern('chlorine', '[Cl]'),
    FragmentPattern('bromine', '[Br]'),
    FragmentPattern('iodine', '[I]'),
    FragmentPattern('lithium', '[Li]'),
    FragmentPattern('sodium', '[Na]'),
    FragmentPattern('potassium', '[K]'),
    FragmentPattern('calcium', '[Ca]'),
    FragmentPattern('magnesium', '[Mg]'),
    FragmentPattern('aluminium', '[Al]'),
    FragmentPattern('barium', '[Ba]'),
    FragmentPattern('bismuth', '[Bi]'),
    FragmentPattern('silver', '[Ag]'),
    FragmentPattern('strontium', '[Sr]'),
    FragmentPattern('zinc', '[Zn]'),
    FragmentPattern('ammonia/ammonium', '[#7]'),
    FragmentPattern('water/hydroxide', '[#8]'),
    FragmentPattern('methyl amine', '[#6]-[#7]'),
    FragmentPattern('sulfide', 'S'),
    FragmentPattern('nitrate', '[#7](=[#8])(-[#8])-[#8]'),
    FragmentPattern('phosphate', '[P](=[#8])(-[#8])(-[#8])-[#8]'),
    FragmentPattern('hexafluorophosphate', '[P](-[#9])(-[#9])(-[#9])(-[#9])(-[#9])-[#9]'),
    FragmentPattern('sulfate', '[S](=[#8])(=[#8])(-[#8])-[#8]'),
    FragmentPattern('methyl sulfonate', '[#6]-[S](=[#8])(=[#8])(-[#8])'),
    FragmentPattern('trifluoromethanesulfonic acid', '[#8]-[S](=[#8])(=[#8])-[#6](-[#9])(-[#9])-[#9]'),
    FragmentPattern('trifluoroacetic acid', '[#9]-[#6](-[#9])(-[#9])-[#6](=[#8])-[#8]'),
    FragmentPattern('1,2-dichloroethane', '[Cl]-[#6]-[#6]-[Cl]'),
    FragmentPattern('1,2-dimethoxyethane', '[#6]-[#8]-[#6]-[#6]-[#8]-[#6]'),
    FragmentPattern('1,4-dioxane', '[#6]-1-[#6]-[#8]-[#6]-[#6]-[#8]-1'),
    FragmentPattern('1-methyl-2-pyrrolidinone', '[#6]-[#7]-1-[#6]-[#6]-[#6]-[#6]-1=[#8]'),
    FragmentPattern('2-butanone', '[#6]-[#6]-[#6](-[#6])=[#8]'),
    FragmentPattern('acetate/acetic acid', '[#8]-[#6](-[#6])=[#8]'),
    FragmentPattern('acetone', '[#6]-[#6](-[#6])=[#8]'),
    FragmentPattern('acetonitrile', '[#6]-[#6]#[N]'),
    FragmentPattern('benzene', '[#6]1[#6][#6][#6][#6][#6]1'),
    FragmentPattern('butanol', '[#8]-[#6]-[#6]-[#6]-[#6]'),
    FragmentPattern('t-butanol', '[#8]-[#6](-[#6])(-[#6])-[#6]'),
    FragmentPattern('chloroform', '[Cl]-[#6](-[Cl])-[Cl]'),
    FragmentPattern('cycloheptane', '[#6]-1-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-1'),
    FragmentPattern('cyclohexane', '[#6]-1-[#6]-[#6]-[#6]-[#6]-[#6]-1'),
    FragmentPattern('dichloromethane', '[Cl]-[#6]-[Cl]'),
    FragmentPattern('diethyl ether', '[#6]-[#6]-[#8]-[#6]-[#6]'),
    FragmentPattern('diisopropyl ether', '[#6]-[#6](-[#6])-[#8]-[#6](-[#6])-[#6]'),
    FragmentPattern('dimethyl formamide', '[#6]-[#7](-[#6])-[#6]=[#8]'),
    FragmentPattern('dimethyl sulfoxide', '[#6]-[S](-[#6])=[#8]'),
    FragmentPattern('ethanol', '[#8]-[#6]-[#6]'),
    FragmentPattern('ethyl acetate', '[#6]-[#6]-[#8]-[#6](-[#6])=[#8]'),
    FragmentPattern('formic acid', '[#8]-[#6]=[#8]'),
    FragmentPattern('heptane', '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]'),
    FragmentPattern('hexane', '[#6]-[#6]-[#6]-[#6]-[#6]-[#6]'),
    FragmentPattern('isopropanol', '[#8]-[#6](-[#6])-[#6]'),
    FragmentPattern('methanol', '[#8]-[#6]'),
    FragmentPattern('N,N-dimethylacetamide', '[#6]-[#7](-[#6])-[#6](-[#6])=[#8]'),
    FragmentPattern('pentane', '[#6]-[#6]-[#6]-[#6]-[#6]'),
    FragmentPattern('propanol', '[#8]-[#6]-[#6]-[#6]'),
    FragmentPattern('pyridine', '[#6]-1=[#6]-[#6]=[#7]-[#6]=[#6]-1'),
    FragmentPattern('t-butyl methyl ether', '[#6]-[#8]-[#6](-[#6])(-[#6])-[#6]'),
    FragmentPattern('tetrahydrofurane', '[#6]-1-[#6]-[#6]-[#8]-[#6]-1'),
    FragmentPattern('toluene', '[#6]-[#6]~1~[#6]~[#6]~[#6]~[#6]~[#6]~1'),
    FragmentPattern('xylene', '[#6]-[#6]~1~[#6](-[#6])~[#6]~[#6]~[#6]~[#6]~1')
)

#: The default value for whether to ensure at least one fragment is left after FragmentRemover is applied.
LEAVE_LAST = True

#: The default value for whether LargestFragmentChooser sees organic fragments as "larger" than inorganic fragments.
PREFER_ORGANIC = False


def is_organic(fragment):
    """Return true if fragment contains at least one carbon atom.

    :param fragment: The fragment as an RDKit Mol object.
    """
    # TODO: Consider a different definition?
    # Could allow only H, C, N, O, S, P, F, Cl, Br, I
    for a in fragment.GetAtoms():
        if a.GetAtomicNum() == 6:
            return True
    return False


class FragmentRemover(object):
    """A class for filtering out fragments using SMARTS patterns."""

    def __init__(self, fragments=REMOVE_FRAGMENTS, leave_last=LEAVE_LAST):
        """Initialize a FragmentRemover with an optional custom list of :class:`~molvs.fragment.FragmentPattern`.

        Setting leave_last to True will ensure at least one fragment is left in the molecule, even if it is matched by a
        :class:`~molvs.fragment.FragmentPattern`. Fragments are removed in the order specified in the list, so place
        those you would prefer to be left towards the end of the list. If all the remaining fragments match the same
        :class:`~molvs.fragment.FragmentPattern`, they will all be left.

        :param fragments: A list of :class:`~molvs.fragment.FragmentPattern` to remove.
        :param bool leave_last: Whether to ensure at least one fragment is left.
        """
        log.debug('Initializing FragmentRemover')
        self.fragments = fragments
        self.leave_last = leave_last

    def __call__(self, mol):
        """Calling a FragmentRemover instance like a function is the same as calling its remove(mol) method."""
        return self.remove(mol)

    def remove(self, mol):
        """Return the molecule with specified fragments removed.

        :param mol: The molecule to remove fragments from.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The molecule with fragments removed.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running FragmentRemover')
        # Iterate FragmentPatterns and remove matching fragments
        for frag in self.fragments:
            # If nothing is left or leave_last and only one fragment, end here
            if mol.GetNumAtoms() == 0 or (self.leave_last and len(Chem.GetMolFrags(mol)) <= 1):
                break
            # Apply removal for this FragmentPattern
            removed = Chem.DeleteSubstructs(mol, frag.smarts, onlyFrags=True)
            if not mol.GetNumAtoms() == removed.GetNumAtoms():
                log.info('Removed fragment: %s', frag.name)
            if self.leave_last and removed.GetNumAtoms() == 0:
                # All the remaining fragments match this pattern - leave them all
                break
            mol = removed
        return mol


class LargestFragmentChooser(object):
    """A class for selecting the largest covalent unit in a molecule with multiple fragments."""

    def __init__(self, prefer_organic=PREFER_ORGANIC):
        """

        If prefer_organic is set to True, any organic fragment will be considered larger than any inorganic fragment. A
        fragment is considered organic if it contains a carbon atom.

        :param bool prefer_organic: Whether to prioritize organic fragments above all others.
        """
        log.debug('Initializing LargestFragmentChooser')
        self.prefer_organic = prefer_organic

    def __call__(self, mol):
        """Calling a LargestFragmentChooser instance like a function is the same as calling its choose(mol) method."""
        return self.choose(mol)

    def choose(self, mol):
        """Return the largest covalent unit.

        The largest fragment is determined by number of atoms (including hydrogens). Ties are broken by taking the
        fragment with the higher molecular weight, and then by taking the first alphabetically by SMILES if needed.

        :param mol: The molecule to choose the largest fragment from.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The largest fragment.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running LargestFragmentChooser')
        # TODO: Alternatively allow a list of fragments to be passed as the mol parameter
        fragments = Chem.GetMolFrags(mol, asMols=True)
        largest = None
        for f in fragments:
            smiles = Chem.MolToSmiles(f, isomericSmiles=True)
            log.debug('Fragment: %s', smiles)
            organic = is_organic(f)
            if self.prefer_organic:
                # Skip this fragment if not organic and we already have an organic fragment as the largest so far
                if largest and largest['organic'] and not organic:
                    continue
                # Reset largest if it wasn't organic and this fragment is organic
                if largest and organic and not largest['organic']:
                    largest = None
            # Count atoms
            atoms = 0
            for a in f.GetAtoms():
                atoms += 1 + a.GetTotalNumHs()
            # Skip this fragment if fewer atoms than the largest
            if largest and atoms < largest['atoms']:
                continue
            # Skip this fragment if equal number of atoms but weight is lower
            weight = rdMolDescriptors.CalcExactMolWt(f)
            if largest and atoms == largest['atoms'] and weight < largest['weight']:
                continue
            # Skip this fragment if equal atoms and equal weight but smiles comes last alphabetically
            if largest and atoms == largest['atoms'] and weight == largest['weight'] and smiles > largest['smiles']:
                continue
            # Otherwise this is the largest so far
            log.debug('New largest fragment: %s (%s)', smiles, atoms)
            largest = {'smiles': smiles, 'fragment': f, 'atoms': atoms, 'weight': weight, 'organic': organic}
        return largest['fragment']



#==============================================================================
#from .normalize import NORMALIZATIONS, MAX_RESTARTS, Normalizer
#==============================================================================
class Normalization(object):
    """A normalization transform defined by reaction SMARTS."""

    def __init__(self, name, transform):
        """
        :param string name: A name for this Normalization
        :param string transform: Reaction SMARTS to define the transformation.
        """
        log.debug('Initializing Normalization: %s', name)
        self.name = name
        self.transform_str = transform

    @memoized_property
    def transform(self):
        log.debug('Loading Normalization transform: %s', self.name)
        return AllChem.ReactionFromSmarts(self.transform_str.encode('utf8'))

    def __repr__(self):
        return 'Normalization({!r}, {!r})'.format(self.name, self.transform_str)

    def __str__(self):
        return self.name
#: The default list of Normalization transforms.        
NORMALIZATIONS = (
    Normalization('Nitro to N+(O-)=O', '[*:1][N,P,As,Sb:2](=[O,S,Se,Te:3])=[O,S,Se,Te:4]>>[*:1][*+1:2]([*-1:3])=[*:4]'),
    Normalization('Sulfone to S(=O)(=O)', '[S+2:1]([O-:2])([O-:3])>>[S+0:1](=[O-0:2])(=[O-0:3])'),
    Normalization('Pyridine oxide to n+O-', '[n:1]=[O:2]>>[n+:1][O-:2]'),
    Normalization('Azide to N=N+=N-', '[*,H:1][N:2]=[N:3]#[N:4]>>[*,H:1][N:2]=[N+:3]=[N-:4]'),
    Normalization('Diazo/azo to =N+=N-', '[*:1]=[N:2]#[N:3]>>[*:1]=[N+:2]=[N-:3]'),
    Normalization('Sulfoxide to -S+(O-)-', '[!O:1][S+0;X3:2](=[O:3])[!O:4]>>[*:1][S+1:2]([O-:3])[*:4]'),
    Normalization('Phosphate to P(O-)=O', '[O,S,Se,Te;-1:1][P+;D4:2][O,S,Se,Te;-1:3]>>[*+0:1]=[P+0;D5:2][*-1:3]'),
    Normalization('Amidinium to C(=NH2+)NH2', '[C,S;X3+1:1]([NX3:2])[NX3!H0:3]>>[*+0:1]([N:2])=[N+:3]'),
    Normalization('Normalize hydrazine-diazonium', '[CX4:1][NX3H:2]-[NX3H:3][CX4:4][NX2+:5]#[NX1:6]>>[CX4:1][NH0:2]=[NH+:3][C:4][N+0:5]=[NH:6]'),
    Normalization('Recombine 1,3-separated charges', '[N,P,As,Sb,O,S,Se,Te;-1:1]-[A:2]=[N,P,As,Sb,O,S,Se,Te;+1:3]>>[*-0:1]=[*:2]-[*+0:3]'),
    Normalization('Recombine 1,3-separated charges', '[n,o,p,s;-1:1]:[a:2]=[N,O,P,S;+1:3]>>[*-0:1]:[*:2]-[*+0:3]'),
    Normalization('Recombine 1,3-separated charges', '[N,O,P,S;-1:1]-[a:2]:[n,o,p,s;+1:3]>>[*-0:1]=[*:2]:[*+0:3]'),
    Normalization('Recombine 1,5-separated charges', '[N,P,As,Sb,O,S,Se,Te;-1:1]-[A+0:2]=[A:3]-[A:4]=[N,P,As,Sb,O,S,Se,Te;+1:5]>>[*-0:1]=[*:2]-[*:3]=[*:4]-[*+0:5]'),
    Normalization('Recombine 1,5-separated charges', '[n,o,p,s;-1:1]:[a:2]:[a:3]:[c:4]=[N,O,P,S;+1:5]>>[*-0:1]:[*:2]:[*:3]:[c:4]-[*+0:5]'),
    Normalization('Recombine 1,5-separated charges', '[N,O,P,S;-1:1]-[c:2]:[a:3]:[a:4]:[n,o,p,s;+1:5]>>[*-0:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
    Normalization('Normalize 1,3 conjugated cation', '[N,O;+0!H0:1]-[A:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]=[*:2]-[*+0:3]'),
    Normalization('Normalize 1,3 conjugated cation', '[n;+0!H0:1]:[c:2]=[N!$(*[O-]),O;+1H0:3]>>[*+1:1]:[*:2]-[*+0:3]'),
    Normalization('Normalize 1,3 conjugated cation', '[N,O;+0!H0:1]-[c:2]:[n!$(*[O-]),o;+1H0:3]>>[*+1:1]=[*:2]:[*+0:3]'),
    Normalization('Normalize 1,5 conjugated cation', '[N,O;+0!H0:1]-[A:2]=[A:3]-[A:4]=[N!$(*[O-]),O;+1H0:5]>>[*+1:1]=[*:2]-[*:3]=[*:4]-[*+0:5]'),
    Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]:[a:2]:[a:3]:[c:4]=[N!$(*[O-]),O;+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]-[*+0:5]'),
    Normalization('Normalize 1,5 conjugated cation', '[N,O;+0!H0:1]-[c:2]:[a:3]:[a:4]:[n!$(*[O-]),o;+1H0:5]>>[*+1:1]=[c:2]:[*:3]:[*:4]:[*+0:5]'),
    Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]1:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]1>>[n+1:1]1:[*:2]:[*:3]:[*:4]:[n+0:5]1'),
    Normalization('Normalize 1,5 conjugated cation', '[n;+0!H0:1]:[a:2]:[a:3]:[a:4]:[n!$(*[O-]);+1H0:5]>>[n+1:1]:[*:2]:[*:3]:[*:4]:[n+0:5]'),
    Normalization('Charge normalization', '[F,Cl,Br,I,At;-1:1]=[O:2]>>[*-0:1][O-:2]'),
    Normalization('Charge recombination', '[N,P,As,Sb;-1:1]=[C+;v3:2]>>[*+0:1]#[C+0:2]'),
)

#: The default value for the maximum number of times to attempt to apply the series of normalizations.
MAX_RESTARTS = 200
# TODO: Rules for canonical resonance / protonation state if a fragment is charged?


class Normalizer(object):
    """A class for applying Normalization transforms.

    This class is typically used to apply a series of Normalization transforms to correct functional groups and
    recombine charges. Each transform is repeatedly applied until no further changes occur.
    """

    def __init__(self, normalizations=NORMALIZATIONS, max_restarts=MAX_RESTARTS):
        """Initialize a Normalizer with an optional custom list of :class:`~molvs.normalize.Normalization` transforms.

        :param normalizations: A list of  :class:`~molvs.normalize.Normalization` transforms to apply.
        :param int max_restarts: The maximum number of times to attempt to apply the series of normalizations (default
                                 200).
        """
        log.debug('Initializing Normalizer')
        self.normalizations = normalizations
        self.max_restarts = max_restarts

    def __call__(self, mol):
        """Calling a Normalizer instance like a function is the same as calling its normalize(mol) method."""
        return self.normalize(mol)

    def normalize(self, mol):
        """Apply a series of Normalization transforms to correct functional groups and recombine charges.

        A series of transforms are applied to the molecule. For each Normalization, the transform is applied repeatedly
        until no further changes occur. If any changes occurred, we go back and start from the first Normalization
        again, in case the changes mean an earlier transform is now applicable. The molecule is returned once the entire
        series of Normalizations cause no further changes or if max_restarts (default 200) is reached.

        :param mol: The molecule to normalize.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The normalized fragment.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running Normalizer')
        # Normalize each fragment separately to get around quirky RunReactants behaviour
        fragments = []
        for fragment in Chem.GetMolFrags(mol, asMols=True):
            fragments.append(self._normalize_fragment(fragment))
        # Join normalized fragments into a single molecule again
        outmol = fragments.pop()
        for fragment in fragments:
            outmol = Chem.CombineMols(outmol, fragment)
        Chem.SanitizeMol(outmol)
        return outmol

    def _normalize_fragment(self, mol):
        for n in xrange(self.max_restarts):
            # Iterate through Normalization transforms and apply each in order
            for normalization in self.normalizations:
                product = self._apply_transform(mol, normalization.transform)
                if product:
                    # If transform changed mol, go back to first rule and apply each again
                    log.info('Rule applied: %s', normalization.name)
                    mol = product
                    break
            else:
                # For loop finishes normally, all applicable transforms have been applied
                return mol
        # If we're still going after max_restarts (default 200), stop and warn, but still return the mol
        log.warn('Gave up normalization after %s restarts', self.max_restarts)
        return mol

    def _apply_transform(self, mol, rule):
        """Repeatedly apply normalization transform to molecule until no changes occur.

        It is possible for multiple products to be produced when a rule is applied. The rule is applied repeatedly to
        each of the products, until no further changes occur or after 20 attempts. If there are multiple unique products
        after the final application, the first product (sorted alphabetically by SMILES) is chosen.
        """
        mols = [mol]
        for n in xrange(20):
            products = {}
            for mol in mols:
                for product in [x[0] for x in rule.RunReactants((mol,))]:
                    if Chem.SanitizeMol(product, catchErrors=True) == 0:
                        products[Chem.MolToSmiles(product, isomericSmiles=True)] = product
            if products:
                mols = [products[s] for s in sorted(products)]
            else:
                # If n == 0, the rule was not applicable and we return None
                return mols[0] if n > 0 else None
#==============================================================================
# from .tautomer import TAUTOMER_TRANSFORMS, TAUTOMER_SCORES, MAX_TAUTOMERS, TautomerCanonicalizer, TautomerEnumerator
#==============================================================================
class TautomerTransform(object):
    """Rules to transform one tautomer to another.

    Each TautomerTransform is defined by a SMARTS pattern where the transform involves moving a hydrogen from the first
    atom in the pattern to the last atom in the pattern. By default, alternating single and double bonds along the
    pattern are swapped accordingly to account for the hydrogen movement. If necessary, the transform can instead define
    custom resulting bond orders and also resulting atom charges.
    """

    BONDMAP = {'-': BondType.SINGLE, '=': BondType.DOUBLE, '#': BondType.TRIPLE, ':': BondType.AROMATIC}
    CHARGEMAP = {'+': 1, '0': 0, '-': -1}

    def __init__(self, name, smarts, bonds=(), charges=()):
        """Initialize a TautomerTransform with a name, SMARTS pattern and optional bonds and charges.

        Specify custom bonds as a string of ``-``, ``=``, ``#``, ``:`` for single, double, triple and aromatic bonds
        respectively. Specify custom charges as ``+``, ``0``, ``-`` for +1, 0 and -1 charges respectively.

        :param string name: A name for this TautomerTransform.
        :param string smarts: SMARTS pattern to match for the transform.
        :param string bonds: Optional specification for the resulting bonds.
        :param string charges: Optional specification for the resulting charges on the atoms.
        """
        self.name = name
        self.tautomer_str = smarts
        self.bonds = [self.BONDMAP[b] for b in bonds]
        self.charges = [self.CHARGEMAP[b] for b in charges]
        # TODO: Raise error (ValueError?) if bonds and charges lists are not the correct length

    @memoized_property
    def tautomer(self):
        return Chem.MolFromSmarts(self.tautomer_str.encode('utf8'))

    def __repr__(self):
        return 'TautomerTransform({!r}, {!r}, {!r}, {!r})'.format(self.name, self.tautomer_str, self.bonds, self.charges)

    def __str__(self):
        return self.name


class TautomerScore(object):
    """A substructure defined by SMARTS and its score contribution to determine the canonical tautomer."""

    def __init__(self, name, smarts, score):
        """Initialize a TautomerScore with a name, SMARTS pattern and score.

        :param name: A name for this TautomerScore.
        :param smarts: SMARTS pattern to match a substructure.
        :param score: The score to assign for this substructure.
        """
        self.name = name
        self.smarts_str = smarts
        self.score = score

    @memoized_property
    def smarts(self):
        return Chem.MolFromSmarts(self.smarts_str.encode('utf8'))

    def __repr__(self):
        return 'TautomerScore({!r}, {!r}, {!r})'.format(self.name, self.smarts_str, self.score)

    def __str__(self):
        return self.name


#: The default list of TautomerTransforms.
TAUTOMER_TRANSFORMS = (
    TautomerTransform('1,3 (thio)keto/enol f', '[CX4!H0][C]=[O,S,Se,Te;X1]'),
    TautomerTransform('1,3 (thio)keto/enol r', '[O,S,Se,Te;X2!H0][C]=[C]'),
    TautomerTransform('1,5 (thio)keto/enol f', '[CX4,NX3;!H0][C]=[C][CH0]=[O,S,Se,Te;X1]'),
    TautomerTransform('1,5 (thio)keto/enol r', '[O,S,Se,Te;X2!H0][CH0]=,:[C][C]=,:[C,N]'),
    TautomerTransform('aliphatic imine f', '[CX4!H0][C]=[NX2]'),
    TautomerTransform('aliphatic imine r', '[NX3!H0][C]=[CX3]'),
    TautomerTransform('special imine f', '[N!H0][C]=[CX3R0]'),
    TautomerTransform('special imine r', '[CX4!H0][c]=,:[n]'),
    TautomerTransform('1,3 aromatic heteroatom H shift f', '[#7!H0][#6R1]=[O,#7X2]'),
    TautomerTransform('1,3 aromatic heteroatom H shift r', '[O,#7;!H0][#6R1]=,:[#7X2]'),
    TautomerTransform('1,3 heteroatom H shift', '[#7,S,O,Se,Te;!H0][#7X2,#6,#15]=[#7,#16,#8,Se,Te]'),
    TautomerTransform('1,5 aromatic heteroatom H shift', '[n,s,o;!H0]:[c,n]:[c]:[c,n]:[n,s,o;H0]'),
    TautomerTransform('1,5 aromatic heteroatom H shift f', '[#7,#16,#8,Se,Te;!H0][#6,nX2]=,:[#6,nX2][#6,#7X2]=,:[#7X2,S,O,Se,Te]'),
    TautomerTransform('1,5 aromatic heteroatom H shift r', '[#7,S,O,Se,Te;!H0][#6,#7X2]=,:[#6,nX2][#6,nX2]=,:[#7,#16,#8,Se,Te]'),
    TautomerTransform('1,7 aromatic heteroatom H shift f', '[#7,#8,#16,Se,Te;!H0][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#6][#6,#7X2]=,:[#7X2,S,O,Se,Te,CX3]'),
    TautomerTransform('1,7 aromatic heteroatom H shift r', '[#7,S,O,Se,Te,CX4;!H0][#6,#7X2]=,:[#6][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[NX2,S,O,Se,Te]'),
    TautomerTransform('1,9 aromatic heteroatom H shift f', '[#7,O;!H0][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#6,#7X2][#6,#7X2]=,:[#7,O]'),
    TautomerTransform('1,11 aromatic heteroatom H shift f', '[#7,O;!H0][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#6,nX2][#6,nX2]=,:[#7X2,O]'),
    TautomerTransform('furanone f', '[O,S,N;!H0][#6X3r5;$([#6][!#6])]=,:[#6X3r5]'),
    TautomerTransform('furanone r', '[#6r5!H0][#6X3r5;$([#6][!#6])]=[O,S,N]'),
    TautomerTransform('keten/ynol f', '[C!H0]=[C]=[O,S,Se,Te;X1]', bonds='#-'),
    TautomerTransform('keten/ynol r', '[O,S,Se,Te;!H0X2][C]#[C]', bonds='=='),
    TautomerTransform('ionic nitro/aci-nitro f', '[C!H0][N+;$([N][O-])]=[O]'),
    TautomerTransform('ionic nitro/aci-nitro r', '[O!H0][N+;$([N][O-])]=[C]'),
    TautomerTransform('oxim/nitroso f', '[O!H0][N]=[C]'),
    TautomerTransform('oxim/nitroso r', '[C!H0][N]=[O]'),
    TautomerTransform('oxim/nitroso via phenol f', '[O!H0][N]=[C][C]=[C][C]=[OH0]'),
    TautomerTransform('oxim/nitroso via phenol r', '[O!H0][c]:[c]:[c]:[c][N]=[OH0]'),
    TautomerTransform('cyano/iso-cyanic acid f', '[O!H0][C]#[N]', bonds='=='),
    TautomerTransform('cyano/iso-cyanic acid r', '[N!H0]=[C]=[O]', bonds='#-'),
    TautomerTransform('formamidinesulfinic acid f', '[O,N;!H0][C][S,Se,Te]=[O]', bonds='=--'),
    TautomerTransform('formamidinesulfinic acid r', '[O!H0][S,Se,Te][C]=[O,N]', bonds='=--'),
    TautomerTransform('isocyanide f', '[C-0!H0]#[N+0]', bonds='#', charges='-+'),
    TautomerTransform('isocyanide r', '[N+!H0]#[C-]', bonds='#', charges='-+'),
    TautomerTransform('phosphonic acid f', '[OH][PH0]', bonds='='),
    TautomerTransform('phosphonic acid r', '[PH]=[O]', bonds='-'),
)

#: The default list of TautomerScores.
TAUTOMER_SCORES = (
    TautomerScore('benzoquinone', '[#6]1([#6]=[#6][#6]([#6]=[#6]1)=,:[N,S,O])=,:[N,S,O]', 25),
    TautomerScore('oxim', '[#6]=[N][OH]', 4),
    TautomerScore('C=O', '[#6]=,:[#8]', 2),
    TautomerScore('N=O', '[#7]=,:[#8]', 2),
    TautomerScore('P=O', '[#15]=,:[#8]', 2),
    TautomerScore('C=hetero', '[#6]=[!#1;!#6]', 1),
    TautomerScore('methyl', '[CX4H3]', 1),
    TautomerScore('guanidine terminal=N', '[#7][#6](=[NR0])[#7H0]', 1),
    TautomerScore('guanidine endocyclic=N', '[#7;R][#6;R]([N])=[#7;R]', 2),
    TautomerScore('aci-nitro', '[#6]=[N+]([O-])[OH]', -4),
)

#: The default value for the maximum number of tautomers to enumerate, a limit to prevent combinatorial explosion.
MAX_TAUTOMERS = 1000


class TautomerCanonicalizer(object):
    """

    """

    def __init__(self, transforms=TAUTOMER_TRANSFORMS, scores=TAUTOMER_SCORES, max_tautomers=MAX_TAUTOMERS):
        """

        :param transforms: A list of TautomerTransforms to use to enumerate tautomers.
        :param scores: A list of TautomerScores to use to choose the canonical tautomer.
        :param max_tautomers: The maximum number of tautomers to enumerate, a limit to prevent combinatorial explosion.
        """
        self.transforms = transforms
        self.scores = scores
        self.max_tautomers = max_tautomers

    def __call__(self, mol):
        """Calling a TautomerCanonicalizer instance like a function is the same as calling its canonicalize(mol) method."""
        return self.canonicalize(mol)

    def canonicalize(self, mol):
        """Return a canonical tautomer by enumerating and scoring all possible tautomers.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The canonical tautomer.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        # TODO: Overload the mol parameter to pass a list of pre-enumerated tautomers
        tautomers = self._enumerate_tautomers(mol)
        if len(tautomers) == 1:
            return tautomers[0]
        # Calculate score for each tautomer
        highest = None
        for t in tautomers:
            smiles = Chem.MolToSmiles(t, isomericSmiles=True)
            log.debug('Tautomer: %s', smiles)
            score = 0
            # Add aromatic ring scores
            ssr = Chem.GetSymmSSSR(t)
            for ring in ssr:
                btypes = {t.GetBondBetweenAtoms(*pair).GetBondType() for pair in pairwise(ring)}
                elements = {t.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring}
                if btypes == {BondType.AROMATIC}:
                    log.debug('Score +100 (aromatic ring)')
                    score += 100
                    if elements == {6}:
                        log.debug('Score +150 (carbocyclic aromatic ring)')
                        score += 150
            # Add SMARTS scores
            for tscore in self.scores:
                for match in t.GetSubstructMatches(tscore.smarts):
                    log.debug('Score %+d (%s)', tscore.score, tscore.name)
                    score += tscore.score
            # Add (P,S,Se,Te)-H scores
            for atom in t.GetAtoms():
                if atom.GetAtomicNum() in {15, 16, 34, 52}:
                    hs = atom.GetTotalNumHs()
                    if hs:
                        log.debug('Score %+d (%s-H bonds)', -hs, atom.GetSymbol())
                        score -= hs
            # Set as highest if score higher or if score equal and smiles comes first alphabetically
            if not highest or highest['score'] < score or (highest['score'] == score and smiles < highest['smiles']):
                log.debug('New highest tautomer: %s (%s)', smiles, score)
                highest = {'smiles': smiles, 'tautomer': t, 'score': score}
        return highest['tautomer']

    @memoized_property
    def _enumerate_tautomers(self):
        return TautomerEnumerator(self.transforms, self.max_tautomers)


class TautomerEnumerator(object):
    """

    """

    def __init__(self, transforms=TAUTOMER_TRANSFORMS, max_tautomers=MAX_TAUTOMERS):
        """

        :param transforms: A list of TautomerTransforms to use to enumerate tautomers.
        :param max_tautomers: The maximum number of tautomers to enumerate (limit to prevent combinatorial explosion).
        """
        self.transforms = transforms
        self.max_tautomers = max_tautomers

    def __call__(self, mol):
        """Calling a TautomerEnumerator instance like a function is the same as calling its enumerate(mol) method."""
        return self.enumerate(mol)

    def enumerate(self, mol):
        """Enumerate all possible tautomers and return them as a list.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: A list of all possible tautomers of the molecule.
        :rtype: list of :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        tautomers = {Chem.MolToSmiles(mol, isomericSmiles=True): copy.deepcopy(mol)}
        done = set()
        while len(tautomers) < self.max_tautomers:
            for tsmiles in sorted(tautomers):
                if tsmiles in done:
                    continue
                for transform in self.transforms:
                    for match in tautomers[tsmiles].GetSubstructMatches(transform.tautomer):
                        # Adjust hydrogens
                        product = copy.deepcopy(tautomers[tsmiles])
                        first = product.GetAtomWithIdx(match[0])
                        last = product.GetAtomWithIdx(match[-1])
                        first.SetNumExplicitHs(max(0, first.GetNumExplicitHs() - 1))
                        last.SetNumExplicitHs(last.GetTotalNumHs() + 1)
                        # Adjust bond orders
                        for bi, pair in enumerate(pairwise(match)):
                            if transform.bonds:
                                product.GetBondBetweenAtoms(*pair).SetBondType(transform.bonds[bi])
                            else:
                                product.GetBondBetweenAtoms(*pair).SetBondType(BondType.DOUBLE if bi % 2 == 0 else BondType.SINGLE)
                        # Adjust charges
                        if transform.charges:
                            for ci, idx in enumerate(match):
                                atom = product.GetAtomWithIdx(idx)
                                atom.SetFormalCharge(atom.GetFormalCharge() + transform.charges[ci])
                        try:
                            Chem.SanitizeMol(product)
                            smiles = Chem.MolToSmiles(product, isomericSmiles=True)
                            log.debug('Applied rule: %s to %s', transform.name, tsmiles)
                            if smiles not in tautomers:
                                log.debug('New tautomer produced: %s' % smiles)
                                tautomers[smiles] = product
                            else:
                                log.debug('Previous tautomer produced again: %s' % smiles)
                        except ValueError:
                            log.debug('ValueError')
                done.add(tsmiles)
            if len(tautomers) == len(done):
                break
        else:
            log.warn('Tautomer enumeration stopped at maximum %s', self.max_tautomers)
        # Clean up stereochemistry
        for tautomer in tautomers.values():
            Chem.AssignStereochemistry(tautomer, force=True, cleanIt=True)
            for bond in tautomer.GetBonds():
                if bond.GetBondType() == BondType.DOUBLE and bond.GetStereo() > BondStereo.STEREOANY:
                    begin = bond.GetBeginAtomIdx()
                    end = bond.GetEndAtomIdx()
                    for othertautomer in tautomers.values():
                        if not othertautomer.GetBondBetweenAtoms(begin, end).GetBondType() == BondType.DOUBLE:
                            neighbours = tautomer.GetAtomWithIdx(begin).GetBonds() + tautomer.GetAtomWithIdx(end).GetBonds()
                            for otherbond in neighbours:
                                if otherbond.GetBondDir() in {BondDir.ENDUPRIGHT, BondDir.ENDDOWNRIGHT}:
                                    otherbond.SetBondDir(BondDir.NONE)
                            Chem.AssignStereochemistry(tautomer, force=True, cleanIt=True)
                            log.debug('Removed stereochemistry from unfixed double bond')
                            break
        return tautomers.values()
#==============================================================================
# from .charge import ACID_BASE_PAIRS, Reionizer, Uncharger
#==============================================================================
class AcidBasePair(object):
    """An acid and its conjugate base, defined by SMARTS.

    A strength-ordered list of AcidBasePairs can be used to ensure the strongest acids in a molecule ionize first.
    """

    def __init__(self, name, acid, base):
        """Initialize an AcidBasePair with the following parameters:

        :param string name: A name for this AcidBasePair.
        :param string acid: SMARTS pattern for the protonated acid.
        :param string base: SMARTS pattern for the conjugate ionized base.
        """
        log.debug('Initializing AcidBasePair: %s', name)
        self.name = name
        self.acid_str = acid
        self.base_str = base

    @memoized_property
    def acid(self):
        log.debug('Loading AcidBasePair acid: %s', self.name)
        return Chem.MolFromSmarts(self.acid_str.encode('utf8'))

    @memoized_property
    def base(self):
        log.debug('Loading AcidBasePair base: %s', self.name)
        return Chem.MolFromSmarts(self.base_str.encode('utf8'))

    def __repr__(self):
        return 'AcidBasePair({!r}, {!r}, {!r})'.format(self.name, self.acid_str, self.base_str)

    def __str__(self):
        return self.name


#: The default list of AcidBasePairs, sorted from strongest to weakest. This list is derived from the Food and Drug
#: Administration Substance Registration System Standard Operating Procedure guide.
ACID_BASE_PAIRS = (
    AcidBasePair('-OSO3H', 'OS(=O)(=O)[OH]', 'OS(=O)(=O)[O-]'),
    AcidBasePair('–SO3H', '[!O]S(=O)(=O)[OH]', '[!O]S(=O)(=O)[O-]'),
    AcidBasePair('-OSO2H', 'O[SD3](=O)[OH]', 'O[SD3](=O)[O-]'),
    AcidBasePair('-SO2H', '[!O][SD3](=O)[OH]', '[!O][SD3](=O)[O-]'),
    AcidBasePair('-OPO3H2', 'OP(=O)([OH])[OH]', 'OP(=O)([OH])[O-]'),
    AcidBasePair('-PO3H2', '[!O]P(=O)([OH])[OH]', '[!O]P(=O)([OH])[O-]'),
    AcidBasePair('-CO2H', 'C(=O)[OH]', 'C(=O)[O-]'),
    AcidBasePair('thiophenol', 'c[SH]', 'c[S-]'),
    AcidBasePair('(-OPO3H)-', 'OP(=O)([OH])[O-]', 'OP(=O)([O-])[O-]'),
    AcidBasePair('(-PO3H)-', '[!O]P(=O)([OH])[O-]', '[!O]P(=O)([O-])[O-]'),
    AcidBasePair('phthalimide', 'O=C2c1ccccc1C(=O)[NH]2', 'O=C2c1ccccc1C(=O)[N-]2'),
    AcidBasePair('CO3H (peracetyl)', 'C(=O)O[OH]', 'C(=O)O[O-]'),
    AcidBasePair('alpha-carbon-hydrogen-nitro group', 'O=N(O)[CH]', 'O=N(O)[C-]'),
    AcidBasePair('-SO2NH2', 'S(=O)(=O)[NH2]', 'S(=O)(=O)[NH-]'),
    AcidBasePair('-OBO2H2', 'OB([OH])[OH]', 'OB([OH])[O-]'),
    AcidBasePair('-BO2H2', '[!O]B([OH])[OH]', '[!O]B([OH])[O-]'),
    AcidBasePair('phenol', 'c[OH]', 'c[O-]'),
    AcidBasePair('SH (aliphatic)', 'C[SH]', 'C[S-]'),
    AcidBasePair('(-OBO2H)-', 'OB([OH])[O-]', 'OB([O-])[O-]'),
    AcidBasePair('(-BO2H)-', '[!O]B([OH])[O-]', '[!O]B([O-])[O-]'),
    AcidBasePair('cyclopentadiene', '[CH2]1C=CC=C1', '[C-]1C=CC=C1'),
    AcidBasePair('-CONH2', 'C(=O)[NH2]', 'C(=O)[NH-]'),
    AcidBasePair('imidazole', 'c1cnc[n]1', 'c1cnc[n-]1'),
    AcidBasePair('-OH', '[CX4][OH]', '[CX4][O-]'),
    AcidBasePair('alpha-carbon-hydrogen-keto group', 'O=C[CH]', 'O=C[C-]'),
    AcidBasePair('alpha-carbon-hydrogen-acetyl ester group', 'OC(=O)[CH]', 'OC(=O)[C-]'),
    AcidBasePair('sp carbon hydrogen', 'C#[CH]', 'C#[C-]'),
    AcidBasePair('alpha-carbon-hydrogen-sulfone group', 'CS(=O)(=O)C[CH]', 'CS(=O)(=O)C[C-]'),
    AcidBasePair('alpha-carbon-hydrogen-sulfoxide group', 'C[SD3](=O)C[CH]', 'C[SD3](=O)C[C-]'),
    AcidBasePair('-NH2', '[CX4][NH2]', '[CX4][NH-]'),
    AcidBasePair('benzyl hydrogen', 'c[CD4H]', 'c[CD3-]'),
    AcidBasePair('sp2-carbon hydrogen', '[CX3]=[CX3H]', '[CX3]=[CX2-]'),
    AcidBasePair('sp3-carbon hydrogen', '[CX4H]', '[CX3-]'),
)


class Reionizer(object):
    """A class to reionize a molecule such that the strongest acids ionize first."""

    def __init__(self, acid_base_pairs=ACID_BASE_PAIRS):
        """Initialize a Reionizer with the following parameter:

        :param acid_base_pairs: A list of :class:`AcidBasePairs <molvs.charge.AcidBasePair>` to reionize, sorted from
                                strongest to weakest.
        """
        log.debug('Initializing Reionizer')
        self.acid_base_pairs = acid_base_pairs

    def __call__(self, mol):
        """Calling a Reionizer instance like a function is the same as calling its reionize(mol) method."""
        return self.reionize(mol)

    def reionize(self, mol):
        """If molecule with multiple acid groups is partially ionized, ensure strongest acids ionize first.

        The algorithm works as follows:

        - Use SMARTS to find the strongest protonated acid and the weakest ionized acid.
        - If the ionized acid is weaker than the protonated acid, swap proton and repeat.

        :param mol: The molecule to reionize.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The reionized molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running Reionizer')
        while True:
            ppos, poccur = self._strongest_protonated(mol)
            ipos, ioccur = self._weakest_ionized(mol)
            if ioccur and poccur and ppos < ipos:
                log.info('Moved proton from %s to %s', self.acid_base_pairs[ppos].name, self.acid_base_pairs[ipos].name)
                patom = mol.GetAtomWithIdx(poccur[-1])
                patom.SetFormalCharge(patom.GetFormalCharge() - 1)
                patom.SetNumExplicitHs(max(0, patom.GetNumExplicitHs() - 1))
                iatom = mol.GetAtomWithIdx(ioccur[-1])
                iatom.SetFormalCharge(iatom.GetFormalCharge() + 1)
                iatom.SetNumExplicitHs(iatom.GetTotalNumHs() + 1)
            else:
                Chem.SanitizeMol(mol)
                return mol

    def _strongest_protonated(self, mol):
        for position, pair in enumerate(self.acid_base_pairs):
            for occurrence in mol.GetSubstructMatches(pair.acid):
                return position, occurrence
        return None, None

    def _weakest_ionized(self, mol):
        for position, pair in enumerate(reversed(self.acid_base_pairs)):
            for occurrence in mol.GetSubstructMatches(pair.base):
                return len(self.acid_base_pairs) - position - 1, occurrence
        return None, None


class Uncharger(object):
    """Class for neutralizing ionized acids and bases.

    This class uncharges molecules by adding and/or removing hydrogens. For zwitterions, hydrogens are moved to
    eliminate charges where possible. However, in cases where there is a positive charge that is not neutralizable, an
    attempt is made to also preserve the corresponding negative charge.

    The method is derived from the neutralise module in `Francis Atkinson's standardiser tool
    <https://github.com/flatkinson/standardiser>`_, which is released under the Apache License v2.0.
    """

    def __init__(self):
        log.debug('Initializing Uncharger')
        #: Neutralizable positive charge (with hydrogens attached)
        self._pos_h = Chem.MolFromSmarts('[+!H0!$(*~[-])]'.encode('utf8'))
        #: Non-neutralizable positive charge (no hydrogens attached)
        self._pos_quat = Chem.MolFromSmarts('[+H0!$(*~[-])]'.encode('utf8'))
        #: Negative charge, not bonded to a positive charge with no hydrogens
        self._neg = Chem.MolFromSmarts('[-!$(*~[+H0])]'.encode('utf8'))
        #: Negative oxygen bonded to [C,P,S]=O, negative aromatic nitrogen?
        self._neg_acid = Chem.MolFromSmarts('[$([O-][C,P,S]=O),$([n-]1nnnc1),$(n1[n-]nnc1)]'.encode('utf8'))

    def __call__(self, mol):
        """Calling an Uncharger instance like a function is the same as calling its uncharge(mol) method."""
        return self.uncharge(mol)

    def uncharge(self, mol):
        """Neutralize molecule by adding/removing hydrogens. Attempts to preserve zwitterions.

        :param mol: The molecule to uncharge.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :return: The uncharged molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        log.debug('Running Uncharger')
        mol = copy.deepcopy(mol)
        # Get atom ids for matches
        p = [x[0] for x in mol.GetSubstructMatches(self._pos_h)]
        q = [x[0] for x in mol.GetSubstructMatches(self._pos_quat)]
        n = [x[0] for x in mol.GetSubstructMatches(self._neg)]
        a = [x[0] for x in mol.GetSubstructMatches(self._neg_acid)]
        # Neutralize negative charges
        if q:
            # Surplus negative charges more than non-neutralizable positive charges
            neg_surplus = len(n) - len(q)
            if a and neg_surplus > 0:
                # zwitterion with more negative charges than quaternary positive centres
                while neg_surplus > 0 and a:
                    # Add hydrogen to first negative acid atom, increase formal charge
                    # Until quaternary positive == negative total or no more negative acid
                    atom = mol.GetAtomWithIdx(a.pop(0))
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                    neg_surplus -= 1
                    log.info('Removed negative charge')
        else:
            #
            for atom in [mol.GetAtomWithIdx(x) for x in n]:
                while atom.GetFormalCharge() < 0:
                    atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
                    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                    log.info('Removed negative charge')
        # Neutralize positive charges
        for atom in [mol.GetAtomWithIdx(x) for x in p]:
            # Remove hydrogen and reduce formal change until neutral or no more hydrogens
            while atom.GetFormalCharge() > 0 and atom.GetNumExplicitHs() > 0:
                atom.SetNumExplicitHs(atom.GetNumExplicitHs() - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                log.info('Removed positive charge')
        return mol

class Standardizer(object):
    """The main class for performing standardization of molecules and deriving parent molecules.

    The primary usage is via the :meth:`~molvs.standardize.Standardizer.standardize` method::

        s = Standardizer()
        mol1 = Chem.MolFromSmiles('C1=CC=CC=C1')
        mol2 = s.standardize(mol1)

    There are separate methods to derive fragment, charge, tautomer, isotope and stereo parent molecules.

    """

    def __init__(self, normalizations=NORMALIZATIONS, acid_base_pairs=ACID_BASE_PAIRS,
                 tautomer_transforms=TAUTOMER_TRANSFORMS, tautomer_scores=TAUTOMER_SCORES,
                 max_restarts=MAX_RESTARTS, max_tautomers=MAX_TAUTOMERS, prefer_organic=PREFER_ORGANIC):
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
        log.debug('Initializing Standardizer')
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
    
    
    def addhs(self,mol):
        from rdkit.Chem import AddHs
        return AddHs(mol)
    
    
    def rmhs(self, mol):
        from rdkit.Chem import RemoveHs
        return RemoveHs(mol)

    def standardize(self, mol):
        """Return a standardized version the given molecule.

        The standardization process consists of the following stages: RDKit
        :rdkit:`RemoveHs <Chem.rdmolops-module.html#RemoveHs>`, RDKit
        :rdkit:`SanitizeMol <Chem.rdmolops-module.html#SanitizeMol>`, :class:`~molvs.metal.MetalDisconnector`,
        :class:`~molvs.normalize.Normalizer`, :class:`~molvs.charge.Reionizer`, RDKit
        :rdkit:`AssignStereochemistry <Chem.rdmolops-module.html#AssignStereochemistry>`.

        :param mol: The molecule to standardize.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :returns: The standardized molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        mol = copy.deepcopy(mol)
        Chem.RemoveHs(mol)
        Chem.SanitizeMol(mol)
        mol = self.disconnect_metals(mol)
        mol = self.normalize(mol)
        mol = self.reionize(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        # TODO: Check this removes symmetric stereocenters
        return mol

    def tautomer_parent(self, mol, skip_standardize=False):
        """Return the tautomer parent of a given molecule.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The tautomer parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        tautomer = self.canonicalize_tautomer(mol)
        tautomer = self.standardize(tautomer)
        return tautomer

    def fragment_parent(self, mol, skip_standardize=False):
        """Return the fragment parent of a given molecule.

        The fragment parent is the largest organic covalent unit in the molecule.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The fragment parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        # TODO: Consider applying FragmentRemover first to remove salts, solvents?
        fragment = self.largest_fragment(mol)
        return fragment

    def stereo_parent(self, mol, skip_standardize=False):
        """Return the stereo parent of a given molecule.

        The stereo parent has all stereochemistry information removed from tetrahedral centers and double bonds.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The stereo parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        else:
            mol = copy.deepcopy(mol)
        mol = Chem.RemoveStereochemistry(mol)
        return mol

    def isotope_parent(self, mol, skip_standardize=False):
        """Return the isotope parent of a given molecule.

        The isotope parent has all atoms replaced with the most abundant isotope for that element.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The isotope parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        else:
            mol = copy.deepcopy(mol)
        # Replace isotopes with common weight
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        return mol

    def charge_parent(self, mol, skip_standardize=False):
        """Return the charge parent of a given molecule.

        The charge parent is the uncharged version of the fragment parent.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The charge parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        # TODO: All ionized acids and bases should be neutralised.
        if not skip_standardize:
            mol = self.standardize(mol)
        fragment = self.fragment_parent(mol, skip_standardize=True)
        if fragment:
            uncharged = self.uncharge(fragment)
            # During final standardization, the Reionizer ensures any remaining charges are in the right places
            uncharged = self.standardize(uncharged)
            return uncharged

    def super_parent(self, mol, skip_standardize=False):
        """Return the super parent of a given molecule.

        THe super parent is fragment, charge, isotope, stereochemistry and tautomer insensitive. From the input
        molecule, the largest fragment is taken. This is uncharged and then isotope and stereochemistry information is
        discarded. Finally, the canonical tautomer is determined and returned.

        :param mol: The input molecule.
        :type mol: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        :param bool skip_standardize: Set to True if mol has already been standardized.
        :returns: The super parent molecule.
        :rtype: :rdkit:`Mol <Chem.rdchem.Mol-class.html>`
        """
        if not skip_standardize:
            mol = self.standardize(mol)
        # We don't need to get fragment parent, because the charge parent is the largest fragment
        mol = self.charge_parent(mol, skip_standardize=True)
        mol = self.isotope_parent(mol, skip_standardize=True)
        mol = self.stereo_parent(mol, skip_standardize=True)
        mol = self.tautomer_parent(mol, skip_standardize=True)
        mol = self.standardize(mol)
        return mol

    def standardize_with_parents(self, mol):
        """"""
        standardized = self.standardize(mol)
        tautomer = self.tautomer_parent(standardized, skip_standardize=True)
        super = self.super_parent(standardized, skip_standardize=True)
        # TODO: Add other parents - have optional argument to specify which are wanted
        mols = {
            'standardized': standardized,
            'tautomer_parent': tautomer,
            'super_parent': super
        }
        return mols

    # TODO: All unique tautomers
    # TODO: All unique fragments (each has to be standardized again?)

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
        return Normalizer(normalizations=self.normalizations, max_restarts=self.max_restarts)

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
    def remove_fragments(self):
        """
        :returns: A callable :class:`~molvs.fragment.FragmentRemover` instance.
        """
        return FragmentRemover()

    @memoized_property
    def largest_fragment(self):
        """
        :returns: A callable :class:`~molvs.fragment.LargestFragmentChooser` instance.
        """
        return LargestFragmentChooser(prefer_organic=self.prefer_organic)

    @memoized_property
    def enumerate_tautomers(self):
        """
        :returns: A callable :class:`~molvs.tautomer.TautomerEnumerator` instance.
        """
        return TautomerEnumerator(transforms=self.tautomer_transforms, max_tautomers=self.max_tautomers)

    @memoized_property
    def canonicalize_tautomer(self):
        """
        :returns: A callable :class:`~molvs.tautomer.TautomerCanonicalizer` instance.
        """
        return TautomerCanonicalizer(transforms=self.tautomer_transforms, scores=self.tautomer_scores,
                                     max_tautomers=self.max_tautomers)

#==============================================================================
# 
#==============================================================================
class StandardizeMol(object):
    """The main class for performing standardization of molecules and deriving parent molecules.

    The primary usage is via the :meth:`~molvs.standardize.Standardizer.standardize` method::

        s = Standardizer()
        mol1 = Chem.MolFromSmiles('C1=CC=CC=C1')
        mol2 = s.standardize(mol1)

    There are separate methods to derive fragment, charge, tautomer, isotope and stereo parent molecules.

    """

    def __init__(self, normalizations=NORMALIZATIONS, acid_base_pairs=ACID_BASE_PAIRS,
                 tautomer_transforms=TAUTOMER_TRANSFORMS, tautomer_scores=TAUTOMER_SCORES,
                 max_restarts=MAX_RESTARTS, max_tautomers=MAX_TAUTOMERS, prefer_organic=PREFER_ORGANIC):
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
        log.debug('Initializing Standardizer')
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
    
    
    def addhs(self,mol):
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
        return Normalizer(normalizations=self.normalizations, max_restarts=self.max_restarts)

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
        return TautomerCanonicalizer(transforms=self.tautomer_transforms, scores=self.tautomer_scores,
                                     max_tautomers=self.max_tautomers)
#==============================================================================
# 
#==============================================================================


def standardize_smiles(smiles):
    """Return a standardized canonical SMILES string given a SMILES string.

    Note: This is a convenience function for quickly standardizing a single SMILES string. It is more efficient to use
    the :class:`~molvs.standardize.Standardizer` class directly when working with many molecules or when custom options
    are needed.

    :param string smiles: The SMILES for the molecule.
    :returns: The SMILES for the standardized molecule.
    :rtype: string.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Standardizer().standardize(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def enumerate_tautomers_smiles(smiles):
    """Return a set of tautomers as SMILES strings, given a SMILES string.

    :param smiles: A SMILES string.
    :returns: A set containing SMILES strings for every possible tautomer.
    :rtype: set of strings.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Standardizer().standardize(mol)
    tautomers = TautomerEnumerator().enumerate(mol)
    return {Chem.MolToSmiles(m, isomericSmiles=True) for m in tautomers}


def canonicalize_tautomer_smiles(smiles):
    """Return a standardized canonical tautomer SMILES string given a SMILES string.

    Note: This is a convenience function for quickly standardizing and finding the canonical tautomer for a single
    SMILES string. It is more efficient to use the :class:`~molvs.standardize.Standardizer` class directly when working
    with many molecules or when custom options are needed.

    :param string smiles: The SMILES for the molecule.
    :returns: The SMILES for the standardize canonical tautomer.
    :rtype: string.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles.encode('utf8'), sanitize=False)
    mol = Standardizer().standardize(mol)
    tautomer = TautomerCanonicalizer().canonicalize(mol)
    return Chem.MolToSmiles(tautomer, isomericSmiles=True)


#==============================================================================
#from .errors import StopValidateError
#==============================================================================
class MolVSError(Exception):
    pass

class StandardizeError(MolVSError):
    pass

class ValidateError(MolVSError):
    pass

class StopValidateError(ValidateError):
    """Called by Validations to stop any further validations from being performed."""
    pass

#==============================================================================
# from .validations import VALIDATIONS
#==============================================================================
class Validation(object):
    """The base class that all :class:`~molvs.validations.Validation` subclasses must inherit from."""

    def __init__(self, log):
        self.log = logging.LoggerAdapter(log, {'validation': type(self).__name__})

    def __call__(self, mol):
        try:
            self.log.debug('Running %s', type(self).__name__)
            self.run(mol)
        except Exception as e:
            if isinstance(e, StopValidateError):
                raise e
            else:
                self.log.debug('Validation failed: %s', e)

    def run(self, mol):
        """"""
        raise NotImplementedError("Validation subclasses must implement the run method")


class SmartsValidation(Validation):
    """Abstract superclass for :class:`Validations <molvs.validations.Validation>` that log a message if a SMARTS
    pattern matches the molecule.

    Subclasses can override the following attributes:
    """

    #: The logging level of the message.
    level = logging.INFO

    #: The message to log if the SMARTS pattern matches the molecule.
    message = 'Molecule matched %(smarts)s'

    #: Whether the SMARTS pattern should match an entire covalent unit.
    entire_fragment = False

    def __init__(self, log):
        super(SmartsValidation, self).__init__(log)
        self._smarts = Chem.MolFromSmarts(self.smarts)

    @property
    def smarts(self):
        """The SMARTS pattern as a string. Subclasses must implement this."""
        raise NotImplementedError('SmartsValidation subclasses must have a smarts attribute')

    def _check_matches(self, mol):
        if mol.HasSubstructMatch(self._smarts):
            self.log.log(self.level, self.message, {'smarts': self.smarts})

    def _check_matches_fragment(self, mol):
        matches = frozenset(frozenset(match) for match in mol.GetSubstructMatches(self._smarts))
        fragments = frozenset(frozenset(frag) for frag in Chem.GetMolFrags(mol))
        if matches & fragments:
            self.log.log(self.level, self.message, {'smarts': self.smarts})

    def run(self, mol):
        if self.entire_fragment:
            self._check_matches_fragment(mol)
        else:
            self._check_matches(mol)


class IsNoneValidation(Validation):
    """Logs an error if ``None`` is passed to the Validator.

    This can happen if RDKit failed to parse an input format. If the molecule is ``None``, no subsequent validations
    will run.
    """

    def run(self, mol):
        if mol is None:
            self.log.error('Molecule is None')
            raise StopValidateError()


class NoAtomValidation(Validation):
    """Logs an error if the molecule has zero atoms.

    If the molecule has no atoms, no subsequent validations will run.
    """

    def run(self, mol):
        if mol.GetNumAtoms() == 0:
            self.log.error('No atoms are present')
            raise StopValidateError()


class DichloroethaneValidation(SmartsValidation):
    """Logs if 1,2-dichloroethane is present.

    This is provided as an example of how to subclass :class:`~molvs.validations.SmartsValidation` to check for the
    presence of a substructure.
    """
    level = logging.INFO
    smarts = '[Cl]-[#6]-[#6]-[Cl]'
    entire_fragment = True
    message = '1,2-Dichloroethane is present'


class FragmentValidation(Validation):
    """Logs if certain fragments are present.

    Subclass and override the ``fragments`` class attribute to customize the list of
    :class:`FragmentPatterns <molvs.fragment.FragmentPattern>`.
    """

    #: A list of :class:`FragmentPatterns <molvs.fragment.FragmentPattern>` to check for.
    fragments = REMOVE_FRAGMENTS

    def run(self, mol):
        for fp in self.fragments:
            matches = frozenset(frozenset(match) for match in mol.GetSubstructMatches(fp.smarts))
            fragments = frozenset(frozenset(frag) for frag in Chem.GetMolFrags(mol))
            if matches & fragments:
                self.log.info('%s is present', fp.name)


class NeutralValidation(Validation):
    """Logs if not an overall neutral system."""

    def run(self, mol):
        charge = Chem.GetFormalCharge(mol)
        if not charge == 0:
            chargestring = '+%s' % charge if charge > 0 else '%s' % charge
            self.log.info('Not an overall neutral system (%s)', chargestring)


class IsotopeValidation(Validation):
    """Logs if molecule contains isotopes."""

    def run(self, mol):
        isotopes = set()
        for atom in mol.GetAtoms():
            isotope = atom.GetIsotope()
            if not isotope == 0:
                isotopes.add('%s%s' % (isotope, atom.GetSymbol()))
        for isotope in isotopes:
            self.log.info('Molecule contains isotope %s', isotope)


#: The default list of :class:`Validations <molvs.validations.Validation>` used by :class:`~molvs.validate.Validator`.
VALIDATIONS = (
    IsNoneValidation,
    NoAtomValidation,
    #DichloroethaneValidation,
    FragmentValidation,
    NeutralValidation,
    IsotopeValidation,
)

#: The default format for log messages.
SIMPLE_FORMAT = '%(levelname)s: [%(validation)s] %(message)s'

#: A more detailed format for log messages. Specify when initializing a Validator.
LONG_FORMAT = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'


class LogHandler(logging.Handler):
    """A simple logging Handler that just stores logs in an array until flushed."""

    def __init__(self):
        logging.Handler.__init__(self)
        self.logs = []

    @property
    def logmessages(self):
        return [self.format(record) for record in self.logs]

    def emit(self, record):
        """Append the record."""
        self.logs.append(record)

    def flush(self):
        """Clear the log records."""
        self.acquire()
        try:
            self.logs = []
        finally:
            self.release()

    def close(self):
        """Close the handler."""
        self.flush()
        logging.Handler.close(self)


class Validator(object):
    """The main class for running :class:`Validations <molvs.validations.Validation>` on molecules."""

    def __init__(self, validations=VALIDATIONS, log_format=SIMPLE_FORMAT, level=logging.INFO, stdout=False, raw=False):
        """Initialize a Validator with the following parameters:

        :param validations: A list of Validations to apply (default: :data:`~molvs.validations.VALIDATIONS`).
        :param string log_format: A string format (default: :data:`~molvs.validate.SIMPLE_FORMAT`).
        :param level: The minimum logging level to output.
        :param bool stdout: Whether to send log messages to standard output.
        :param bool raw: Whether to return raw :class:`~logging.LogRecord` objects instead of formatted log strings.
        """
        self.raw = raw
        # Set up logger and add default LogHandler
        self.log = logging.getLogger(type(self).__name__)
        self.log.setLevel(level)
        self.handler = LogHandler()
        self.handler.setFormatter(logging.Formatter(log_format))
        self.log.addHandler(self.handler)
        # Add stdout StreamHandler if specified in parameters
        if stdout:
            strhdlr = logging.StreamHandler(sys.stdout)
            strhdlr.setFormatter(logging.Formatter(log_format))
            self.log.addHandler(strhdlr)
        # Instantiate the validations
        self.validations = [validation(self.log) for validation in validations]

    def __call__(self, mol):
        """Calling a Validator instance like a function is the same as calling its
        :meth:`~molvs.validate.Validator.validate` method."""
        return self.validate(mol)

    def validate(self, mol):
        """"""
        # Clear any log messages from previous runs
        self.handler.flush()
        # Run every validation, stopping if StopValidateError is raised
        for validation in self.validations:
            try:
                validation(mol)
            except StopValidateError:
                break
        return self.handler.logs if self.raw else self.handler.logmessages


def validate_smiles(smiles):
    """Return log messages for a given SMILES string using the default validations.

    Note: This is a convenience function for quickly validating a single SMILES string. It is more efficient to use
    the :class:`~molvs.validate.Validator` class directly when working with many molecules or when custom options
    are needed.

    :param string smiles: The SMILES for the molecule.
    :returns: A list of log messages.
    :rtype: list of strings.
    """
    # Skip sanitize as standardize does this anyway
    mol = Chem.MolFromSmiles(smiles.encode('utf8'))
    logs = Validator().validate(mol)
    return logs


if __name__ == '__main__':
    mol = Chem.MolFromSmiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
    from rdkit.Chem import AddHs
    
    mol = AddHs(mol)
    s = Standardizer()
    smol = s.standardize(mol)
    standardize_smiles('[Na]OC(=O)c1ccc(C[S+2]([O-])([O-]))cc1')
    l = LargestFragmentChooser()
    lmol = l.choose(mol)
    print(Chem.MolToSmiles(lmol, isomericSmiles=True))
    print (standardize_smiles('C[n+]1c([N-](C))cccc1'))











