![](https://img.shields.io/pypi/l/Django.svg) ![](https://img.shields.io/badge/dependencies-rdkit%2C%20pybel-green.svg) ![](https://img.shields.io/badge/platforms-linux%2C%20windows-brightgreen.svg)<br>
## Introduction
PyBioMed provides various user-friendly and highly customized APIs to calculate various features(descriptors) of biological molecules(chemicals, proteins and DNA/RNA descriptors) and complex interaction samples conveniently, which aims at building integrated analysis pipelines from data acquisition, data checking, and descriptor calculation to modeling.
## Installation
### Install Pybel and RDKit
To install Pybel: http://openbabel.org/docs/current/UseTheLibrary/PythonInstall.html<br>
To install RDKit: http://www.rdkit.org/docs/Install.html<br>
### Install PyBioMed
PyBioMed has been successfully tested on Linux and Windows systems. After installing RDKit and pybel successfully, The author could download the PyBioMed package via: https://github.com/gadsbyfly/PyBioMed/blob/master/PyBioMed/download/PyBioMed-1.0.zip. <br>
The install process of PyBioMed is very easy:

#### On Windows:

(1): download the PyBioMed-1.0.zip

(2): extract or uncompress the PyBioMed-1.0.zip file

(3): cd PyBioMed-1.0

(4): python setup.py install

#### On Linux:

(1): download the PyBioMed package (.zip)

(2): unzip PyBioMed-1.0.zip

(3): cd PyBioMed-1.0

(4): python setup.py install or sudo python setup.py install
### Recommended installation
There is a simplest configuration on Ubuntu 14.04(Just click mouse to install): <br>
(1) use the ubuntu software center to search for 'Synaptic Package Manager' and install it. Use 'Synaptic Package Manager' to search for 'openbabel, libopenbabel4, python-openbabel, libopenbabel-dev' and then install them. Search for 'python-rdkit, librdkit1, rdkit-data' and install them. This will make sure the right installation of Pybel and RDKit.<br>
(2) Download PyBioMed-1.0.zip, 'cd PyBioMed-1.0', 'python setup.py install or sudo python setup.py install'.
## Documentation
(1)The online version of the documentation is available here:<br>
http://projects.scbdd.com/pybiomed/index.html <br>
or http://pybiomed.readthedocs.io/en/latest/<br>
(2)Quick start examples: http://projects.scbdd.com/pybiomed/User_guide.html<br>
or http://pybiomed.readthedocs.io/en/latest/User_guide.html<br>
(3)Application examples(pipelines): http://projects.scbdd.com/pybiomed/application.html<br>
or http://pybiomed.readthedocs.io/en/latest/application.html
## Contact
If you have questions or suggestions, please contact:
gadsby@163.com, biomed@csu.edu.cn and oriental-cds@163.com

Please see the file LICENSE for details about the "New BSD"
license which covers this software and its associated data and
documents.

# Copyright (C) 2015-2017 CBDD Group
