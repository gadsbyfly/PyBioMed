"""This file was a duplicate of PyBioMed.PyPretreat.PyPretreatDNA."""

# Core Library modules
import warnings

# First party modules
from PyBioMed.PyPretreat.PyPretreatDNA import *  # noqa

warnings.warn(
    "The PyBioMed.PyPretreat.PyDNAutil is deprecated "
    "and will be removed in 2021. "
    "Use PyBioMed.PyPretreat.PyPretreatDNA instead.",
    DeprecationWarning,
)
