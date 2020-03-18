# Core Library modules
import sys

# Third party modules
from setuptools import setup

packagedata = {
    "PyBioMed": [
        "PyDNA/*",
        "PyGetMol/*",
        "PyInteraction/*",
        "PyPretreat/*",
        "PyProtein/*",
        "PyMolecule/*",
        "PyMolecule/*",
        "test/*",
        "test/test_data/*",
        "example/caco2/*",
        "example/dna/*",
        "example/dpi/*",
        "example/solubility/*",
        "example/subcell/*",
        "doc/html/*",
        "doc/tml/_sources/*",
        "doc/html/_sources/reference/*",
        "doc/html/_images/*",
        "doc/html/_modules/*",
        "doc/html/_static/*",
        "doc/html/reference/*",
        "doc/*",
        "doc/Descriptor/*",
    ]
}

setup(package_data=packagedata, package_dir={"PyBioMed": "PyBioMed"})
