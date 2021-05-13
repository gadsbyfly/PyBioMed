# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################
The calculation of molecular topological indices based on its topological

structure. You can get 25 molecular topological descriptors. You can freely

use and distribute it. If you hava  any problem, you could contact with us timely!

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

# Standard library
from collections import Counter
# Third party modules
import numpy
import scipy
from rdkit import Chem
from rdkit.Chem import GraphDescriptors as GD
from rdkit.Chem import rdchem

periodicTable = rdchem.GetPeriodicTable()

Version = 1.0
################################################################


def _GetPrincipleQuantumNumber(atNum):
    """
    #################################################################
    *Internal Use Only*

    Get the principle quantum number of atom with atomic

    number equal to atNum
    #################################################################
    """
    if atNum <= 2:
        return 1
    elif atNum <= 10:
        return 2
    elif atNum <= 18:
        return 3
    elif atNum <= 36:
        return 4
    elif atNum <= 54:
        return 5
    elif atNum <= 86:
        return 6
    else:
        return 7


def CalculateWeiner(mol):
    """
    #################################################################
    Calculation of Weiner number in a molecule

    ---->W

    Usage:

        result=CalculateWeiner(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    return 1.0 / 2 * sum(sum(Chem.GetDistanceMatrix(mol)))


def CalculateMeanWeiner(mol):

    """
    #################################################################
    Calculation of Mean Weiner number in a molecule

    ---->AW

    Usage:

        result=CalculateWeiner(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    N = mol.GetNumAtoms()
    WeinerNumber = CalculateWeiner(mol)
    return 2.0 * WeinerNumber / (N * (N - 1))


def CalculateBalaban(mol):

    """
    #################################################################
    Calculation of Balaban index in a molecule

    ---->J

    Usage:

        result=CalculateBalaban(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    adjMat = Chem.GetAdjacencyMatrix(mol)
    Distance = Chem.GetDistanceMatrix(mol)
    Nbond = mol.GetNumBonds()
    Natom = mol.GetNumAtoms()
    S = numpy.sum(Distance, axis=1)
    mu = Nbond - Natom + 1
    sumk = 0.0
    for i in range(len(Distance)):
        si = S[i]
        for j in range(i, len(Distance)):
            if adjMat[i, j] == 1:
                sumk += 1.0 / numpy.sqrt(si * S[j])
    if mu + 1 != 0:
        J = float(Nbond) / float(mu + 1) * sumk
    else:
        J = 0
    return J


def CalculateGraphDistance(mol):
    """
    #################################################################
    Calculation of graph distance index

    ---->Tigdi(log value)

    Usage:

        result=CalculateGraphDistance(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    Distance= Chem.GetDistanceMatrix(mol) 
    temp =  Counter(Distance.flatten())
    res = sum(
        map(lambda x: (x*0.5)**2, temp.values())
    )
    return numpy.log10(res)


def CalculateDiameter(mol):
    """
    #################################################################
    Calculation of diameter, which is 	Largest value

    in the distance matrix [Petitjean 1992].

    ---->diametert

    Usage:

        result=CalculateDiameter(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    Distance = Chem.GetDistanceMatrix(mol)

    return Distance.max()


def CalculateRadius(mol):
    """
    #################################################################
    Calculation of radius based on topology.

    It is :If ri is the largest matrix entry in row i of the distance

    matrix D,then the radius is defined as the smallest of the ri

    [Petitjean 1992].

    ---->radiust

    Usage:

        result=CalculateRadius(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    Distance = Chem.GetDistanceMatrix(mol)
    temp = []
    for i in Distance:
        temp.append(max(i))
    return min(temp)


def CalculatePetitjean(mol):
    """
    #################################################################
    Calculation of Petitjean based on topology.

    Value of (diameter - radius) / diameter as defined in [Petitjean 1992].

    ---->petitjeant

    Usage:

        result=CalculatePetitjean(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    diameter = CalculateDiameter(mol)
    radius = CalculateRadius(mol)
    return 1 - radius / float(diameter)


def CalculateXuIndex(mol):
    """
    #################################################################
    Calculation of Xu index

    ---->Xu

    Usage:

        result=CalculateXuIndex(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    sigma = scipy.sum(Distance, axis=1)
    temp1 = 0.0
    temp2 = 0.0
    for i in range(nAT):
        temp1 = temp1 + deltas[i] * ((sigma[i]) ** 2)
        temp2 = temp2 + deltas[i] * (sigma[i])
    Xu = numpy.sqrt(nAT) * numpy.log(temp1 / temp2)

    return Xu


def CalculateGutmanTopo(mol):
    """
    #################################################################
    Calculation of Gutman molecular topological index based on

    simple vertex degree

    ---->GMTI(log value)

    Usage:

        result=CalculateGutmanTopo(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    nAT = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    Distance = Chem.GetDistanceMatrix(mol)
    res = 0.0
    for i in range(nAT):
        for j in range(i + 1, nAT):
            res = res + deltas[i] * deltas[j] * Distance[i, j]

    return numpy.log10(res)

def _HKDeltas(mol,skipHs=1):
    """
    #################################################################
    *Internal Use Only*
    
    Calculation of modified delta value for a molecule
    
    res---->list type
    #################################################################
    """
    global periodicTable
    res=[]
    for atom in mol.GetAtoms():
        n=atom.GetAtomicNum()
        if n>1:
            nV=periodicTable.GetNOuterElecs(n)
            nHs=atom.GetTotalNumHs()
            if n<10:
                res.append(float(nV-nHs))
            else:
                res.append(float(nV-nHs)/float(n-nV-1))
        elif not skipHs:
            res.append(0.0)
    return res

def CalculateGutmanVTopo(mol):
    """
    #################################################################
    Calculation of Gutman molecular topological index based on
    
    valence vertex degree(log1o)
    
    ---->GMTIV
    
    Usage: 
        
        result=CalculateGutmanVTopo(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    nAT=mol.GetNumAtoms()
    deltas=_HKDeltas(mol)
    Distance= Chem.GetDistanceMatrix(mol)
    res=0.0
    for i in range(nAT):
        for j in range(i+1,nAT):
            res=res+deltas[i]*deltas[j]*Distance[i,j]

    return numpy.log10(res)


def CalculatePolarityNumber(mol):
    """
    #################################################################
    Calculation of Polarity number.

    It is the number of pairs of vertexes at

    distance matrix equal to 3

    ---->Pol

    Usage:

        result=CalculatePolarityNumber(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    Distance = Chem.GetDistanceMatrix(mol)
    res = 1.0 / 2 * sum(sum(Distance == 3))

    return res


def CalculatePoglianiIndex(mol):
    """
    #################################################################
    Calculation of Poglicani index

    The Pogliani index (Dz) is the sum over all non-hydrogen atoms

    of a modified vertex degree calculated as the ratio

    of the number of valence electrons over the principal

    quantum number of an atom [L. Pogliani, J.Phys.Chem.

    1996, 100, 18065-18077].

    ---->DZ

    Usage:

        result=CalculatePoglianiIndex(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    res = 0.0
    for atom in mol.GetAtoms():
        n = atom.GetAtomicNum()
        nV = periodicTable.GetNOuterElecs(n)
        mP = _GetPrincipleQuantumNumber(n)
        res = res + (nV + 0.0) / mP
    return res


def CalculateIpc(mol):

    """
    #################################################################
    This returns the information content of the coefficients of the

    characteristic polynomial of the adjacency matrix of a

    hydrogen-suppressed graph of a molecule.

    'avg = 1' returns the information content divided by the total

    population.

    From D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67,

    4517-4533 (1977)

     ---->Ipc(log value)

    Usage:

        result=CalculateIpc(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    temp = GD.Ipc(mol)
    if temp > 0:
        return numpy.log10(temp)
    else:
        return "NaN"


def CalculateBertzCT(mol):
    """
    #################################################################
    A topological index meant to quantify "complexity" of molecules.

    Consists of a sum of two terms, one representing the complexity

    of the bonding, the other representing the complexity of the

    distribution of heteroatoms.

    From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)

    ---->BertzCT(log value)

    Usage:

        result=CalculateBertzCT(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    temp = GD.BertzCT(mol)
    if temp > 0:
        return numpy.log10(temp)
    else:
        return "NaN"


def CalculateHarary(mol):
    """
    #################################################################
    Calculation of Harary number

    ---->Thara

    Usage:

        result=CalculateHarary(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """

    Distance = numpy.array(Chem.GetDistanceMatrix(mol), "d")

    return 1.0 / 2 * (sum(1.0 / Distance[Distance != 0]))


def CalculateSchiultz(mol):
    """
    #################################################################
    Calculation of Schiultz number

    ---->Tsch(log value)

    Usage:

        result=CalculateSchiultz(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    Distance = numpy.array(Chem.GetDistanceMatrix(mol), "d")
    Adjacent = numpy.array(Chem.GetAdjacencyMatrix(mol), "d")
    VertexDegree = sum(Adjacent)

    return sum(scipy.dot((Distance + Adjacent), VertexDegree))


def CalculateZagreb1(mol):
    """
    #################################################################
    Calculation of Zagreb index with order 1 in a molecule

    ---->ZM1

    Usage:

        result=CalculateZagreb1(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """

    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    return sum(numpy.array(deltas) ** 2)


def CalculateZagreb2(mol):

    """
    #################################################################
    Calculation of Zagreb index with order 2 in a molecule

    ---->ZM2

    Usage:

        result=CalculateZagreb2(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    ke = [
        x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree()
        for x in mol.GetBonds()
    ]

    return sum(ke)


def CalculateMZagreb1(mol):
    """
    #################################################################
    Calculation of Modified Zagreb index with order 1 in a molecule

    ---->MZM1

    Usage:

        result=CalculateMZagreb1(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")
    res = sum((1.0 / deltas) ** 2)
    return res


def CalculateMZagreb2(mol):
    """
    #################################################################
    Calculation of Modified Zagreb index with order 2 in a molecule

    ---->MZM2

    Usage:

        result=CalculateMZagreb2(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    cc = [
        x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree()
        for x in mol.GetBonds()
    ]
    while 0 in cc:
        cc.remove(0)
    cc = numpy.array(cc, "d")
    res = sum((1.0 / cc) ** 2)
    return res


def CalculateQuadratic(mol):
    """
    #################################################################
    Calculation of Quadratic index in a molecule

    ---->Qindex

    Usage:

        result=CalculateQuadratic(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    M = CalculateZagreb1(mol)
    N = mol.GetNumAtoms()
    return 3 - 2 * N + M / 2.0


def CalculatePlatt(mol):
    """
    #################################################################
    Calculation of Platt number in a molecule

    ---->Platt

    Usage:

        result=CalculatePlatt(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    cc = [
        x.GetBeginAtom().GetDegree() + x.GetEndAtom().GetDegree() - 2
        for x in mol.GetBonds()
    ]
    return sum(cc)


def CalculateSimpleTopoIndex(mol):
    """
    #################################################################
    Calculation of the logarithm of the simple topological index by Narumi,

    which is defined as the product of the vertex degree.

    ---->Sito

    Usage:

        result=CalculateSimpleTopoIndex(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")

    res = numpy.prod(deltas)
    if res > 0:
        return numpy.log10(res)
    else:
        return "NaN"


def CalculateHarmonicTopoIndex(mol):
    """
    #################################################################
    Calculation of harmonic topological index proposed by Narnumi.

    ---->Hato

    Usage:

        result=CalculateHarmonicTopoIndex(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")
    nAtoms = mol.GetNumAtoms()

    res = nAtoms / sum(1.0 / deltas)

    return res


def CalculateGeometricTopoIndex(mol):
    """
    #################################################################
    Geometric topological index by Narumi

    ---->Geto

    Usage:

        result=CalculateGeometricTopoIndex(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    nAtoms = mol.GetNumAtoms()
    deltas = [x.GetDegree() for x in mol.GetAtoms()]
    while 0 in deltas:
        deltas.remove(0)
    deltas = numpy.array(deltas, "d")

    temp = numpy.prod(deltas)
    res = numpy.power(temp, 1.0 / nAtoms)

    return res


def CalculateArithmeticTopoIndex(mol):
    """
    #################################################################
    Arithmetic topological index by Narumi

    ---->Arto

    Usage:

        result=CalculateArithmeticTopoIndex(mol)

        Input: mol is a molecule object

        Output: result is a numeric value
    #################################################################
    """
    nAtoms = mol.GetNumAtoms()
    nBonds = mol.GetNumBonds()

    res = 2.0 * nBonds / nAtoms
    return res

def CalculateDistanceEqualityTotalInf(mol):
    
    """
    #################################################################
    Total information index on distance equality
    
    -->DET
    
    Usage: 
        
        result=CalculateDistanceEqualityTotalInf(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    Distance= Chem.GetDistanceMatrix(mol)
    nAT=mol.GetNumAtoms()
    n=1./2*nAT**2-nAT
    DisType=int(Distance.max())
    res=0.0
    #res=numpy.zeros(DisType,numpy.float)
    for i in range(DisType):
        cc=1./2*sum(sum(Distance==i+1))
        res += cc*numpy.log2(cc)

    return n*numpy.log2(n)-res


def CalculateDistanceEqualityMeanInf(mol):
    
    """
    #################################################################
    Mean information index on distance equality
    
    -->IDE
    
    Usage: 
        
        result=CalculateDistanceEqualityMeanInf(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    Distance= Chem.GetDistanceMatrix(mol)

    nAT=mol.GetNumAtoms()
    n=1./2*nAT**2-nAT
    DisType=int(Distance.max())
    res=0.0
    cc=numpy.zeros(DisType,numpy.float)
    for i in range(DisType):
        cc[i]=1./2*sum(sum(Distance==i+1))

    res=_CalculateEntropy(cc/n)
      
    return res

def CalculateMolSizeTotalInf(mol):
    """
    #################################################################
    Total information index on molecular size
    
    -->ISIZ
    
    Usage: 
        
        result=CalculateMolSizeTotalInf(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    nAT=Hmol.GetNumAtoms()
    ISIZ=nAT*numpy.log2(nAT)
    return ISIZ

def _CalculateEntropy(Probability):
    """
    #################################################################
    **Internal used only**
    
    Calculation of entropy (Information content) for probability given
    #################################################################
    """
    res=0.0
    for i in Probability:
        if i!=0:       
            res=res-i*numpy.log2(i)

    return res

def CalculateVertexEqualityTotalInf(mol):
    """
    #################################################################
    Total information index on vertex equality
    
    -->IVDE
    
    Usage: 
        
        result=CalculateVertexEqualityTotalInf(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    deltas=[x.GetDegree() for x in mol.GetAtoms()]
    res=0.0
    while 0 in deltas:
        deltas.remove()
    for i in range(max(deltas)):
        cc=deltas.count(i+1)
        if cc==0:
            res=res
        else:
            res += cc*numpy.log2(cc)

    n=len(deltas)

    return n*numpy.log2(n)-res

def CalculateAtomCompTotalInf(mol):
    """
    #################################################################
    Total information index on atomic composition
    
    -->TIAC
    
    Usage: 
        
        result=CalculateAtomCompTotalInf(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    Hmol=Chem.AddHs(mol)
    nAtoms=Hmol.GetNumAtoms()
    IC=[]
    for i in range(nAtoms):
        at=Hmol.GetAtomWithIdx(i)
        IC.append(at.GetAtomicNum())
    Unique=numpy.unique(IC)
    NAtomType=len(Unique)
    #NTAtomType=numpy.zeros(NAtomType,numpy.float)
    res=0.0
    for i in range(NAtomType):
        cc=IC.count(Unique[i])
        res += cc*numpy.log2(cc)

    if nAtoms!=0:
        return nAtoms*numpy.log2(nAtoms)-res
    else:
        return 0.0

    
def CalculateGravitationalTopoIndex(mol):
    """
    #################################################################
    Gravitational topological index based on topological distance 
    
    instead of intermolecular distance.
    
    ---->Gravto
    
    Usage: 
        
        result=CalculateGravitationalTopoIndex(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    nAT=mol.GetNumAtoms()
    Distance= Chem.GetDistanceMatrix(mol) 
    res=0.0
    Atom=mol.GetAtoms()
    for i in range(nAT-1):
        for j in range(i+1,nAT):
            temp=Atom[i].GetMass()*Atom[j].GetMass()
            res=res+temp/numpy.power(Distance[i][j],2)
    
    return res/100      

def CalculateSimpleTopovIndex(mol):
    """
    #################################################################
    Calculation of the logarithm of the simple topological index by Narumi,
    
    which is defined as the product of the vertex degree.
    
    ---->Sitov
    
    Usage: 
        
        result=CalculateSimpleTopovIndex(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    deltas=_HKDeltas(mol,skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    
    res=numpy.prod(deltas)
    
    return numpy.log(res)

def CalculateGeometricTopovIndex(mol):
    """
    #################################################################
    Geometric topological index by Narumi
    
    ---->Getov
    
    Usage: 
        
        result=CalculateGeometricTopovIndex(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    nAtoms=mol.GetNumAtoms()
    deltas=_HKDeltas(mol,skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')
    
    temp=numpy.prod(deltas)
    res=numpy.power(temp,1./nAtoms)

    return res 

def CalculateHarmonicTopovIndex(mol):
    """
    #################################################################
    Calculation of harmonic topological index proposed by Narnumi.
    
    ---->Hatov
    
    Usage: 
        
        result=CalculateHarmonicTopovIndex(mol)
        
        Input: mol is a molecule object
        
        Output: result is a numeric value
    #################################################################
    """
    deltas=_HKDeltas(mol,skipHs=0)
    while 0 in deltas:
        deltas.remove(0)
    deltas=numpy.array(deltas,'d')  
    nAtoms=mol.GetNumAtoms()
    
    res=nAtoms/sum(1./deltas)
    
    return res



_Topology = {
    "W": CalculateWeiner,
    "AW": CalculateMeanWeiner,
    "J": CalculateBalaban,
    "Tigdi": CalculateGraphDistance,
    "Xu": CalculateXuIndex,
    "GMTI": CalculateGutmanTopo,
    "Pol": CalculatePolarityNumber,
    "DZ": CalculatePoglianiIndex,
    "Ipc": CalculateIpc,
    "BertzCT": CalculateBertzCT,
    "Thara": CalculateHarary,
    "Tsch": CalculateSchiultz,
    "ZM1": CalculateZagreb1,
    "ZM2": CalculateZagreb2,
    "MZM1": CalculateMZagreb1,
    "MZM2": CalculateMZagreb2,
    "Qindex": CalculateQuadratic,
    "Platt": CalculatePlatt,
    "diametert": CalculateDiameter,
    "radiust": CalculateRadius,
    "petitjeant": CalculatePetitjean,
    "Sito": CalculateSimpleTopoIndex,
    "Hato": CalculateHarmonicTopoIndex,
    "Geto": CalculateGeometricTopoIndex,
    "Arto": CalculateArithmeticTopoIndex,
    'GMTIV':CalculateGutmanVTopo,
    'IDET':CalculateDistanceEqualityTotalInf,
    'IDE':CalculateDistanceEqualityMeanInf,
    'ISIZ':CalculateMolSizeTotalInf,
    'IVDE':CalculateVertexEqualityTotalInf,
    'TIAC':CalculateAtomCompTotalInf,
    'Gravto':CalculateGravitationalTopoIndex,
    'Hatov':CalculateHarmonicTopovIndex,
    'Sitov':CalculateSimpleTopovIndex,
    'Getov':CalculateGeometricTopovIndex,
}


def GetTopology(mol):
    """
    #################################################################
    Get the dictionary of constitutional descriptors for given

    moelcule mol

    Usage:

        result=CalculateWeiner(mol)

        Input: mol is a molecule object

        Output: result is a dict form containing all topological indices.
    #################################################################
    """
    result = {}
    for DesLabel in _Topology.keys():
        result[DesLabel] = round(_Topology[DesLabel](mol), 3)
    return result


def _GetHTMLDoc():
    """
    #################################################################
    Write HTML documentation for this module.
    #################################################################
    """
    import pydoc

    pydoc.writedoc("topology")


################################################################################
#####################################################################

if __name__ == "__main__":

    smis = ["CCCC", "CCCCC", "CCCCCC", "CC(N)C(=O)O", "CC(N)C(=O)[O-]"]
    for index, smi in enumerate(smis):
        m = Chem.MolFromSmiles(smi)
        print(index + 1)
        print(smi)
        print("\t", GetTopology(m))
        print("\t", len(GetTopology(m)))
