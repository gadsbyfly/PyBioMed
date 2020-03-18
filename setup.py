# Core Library modules
import sys
from distutils.core import setup

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

setup(
    name="PyBioMed",
    version="1.0",
    description="biological molecular representation Python package (PyBioMed)",
    long_description="PyBioMed is a feature-rich packages used for the characterization of various complex biological molecules and interaction samples, such as chemicals, proteins, DNA, and their interactions. PyBioMed calculates nine types of features including chemical descriptors or molecular fingerprints, structural and physicochemical features of proteins and peptides from amino acid sequence, composition and physicochemical features of DNA from their primary sequences, chemical-chemical interaction features, chemical-protein interaction features, chemical-DNA interaction features, protein-protein interaction features, protein-DNA interaction features, and DNA-DNA interaction features. PyBioMed can also pretreating molecule structures, protein sequences and DNA sequence. In order to be convenient to users, PyBioMed provides the module to get molecule structures, protein sequence and DNA sequence from Internet.",
    author="CBDD Group",
    author_email="gadsby@163.com and oriental-cds@163.com",
    license="BSD",
    packages=["PyBioMed"],
    zip_safe=False,
    package_data=packagedata,
    # data_files = datafiles,
    package_dir={"PyBioMed": "PyBioMed"},
    scripts=[],
    py_modules=[],
)
