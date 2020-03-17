*******
Testing
*******

Requirements for testing
========================
PyBioMed requires RDKit and pybel packages.
If you don't already have the packages installed, follow
the directions here
https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html

http://www.rdkit.org/

Testing an installed package
============================

If you have a file-based (not a Python egg) installation you can
test the installed package with

>>> from PyBioMed.test import test_PyBioMed
>>> test_PyBioMed.test_pybiomed()
