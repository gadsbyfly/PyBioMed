# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
##############################################################################

A class used for computing different types of DNA descriptors! 

You can freely use and distribute it. If you have any problem, 

you could contact with us timely.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.07.11

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""

from PyDNAutil import GetData, GeneratePhycheValue
from functools import reduce


def CheckAcc(lag, k):
    """
    #################################################################
    Check ACC parameter validation.
    #################################################################
    """
    try:
        if not isinstance(lag, int) or lag <= 0:
            raise ValueError("Error, parameter lag must be an int type and larger than 0.")
        elif not isinstance(k, int) or lag <= 0:
            raise ValueError("Error, parameter k must be an int type and larger than 0.")
    except ValueError:
        raise


def ReadyAcc(input_data, k, phyche_index=None, all_property=False, extra_phyche_index=None):
    """
    #################################################################
    Public function for get sequence_list and phyche_value.
    #################################################################
    """
    sequence_list = GetData(input_data)
    if phyche_index is None:
        phyche_index = []
    if extra_phyche_index is None:
        extra_phyche_index = {}
    phyche_value = GeneratePhycheValue(k, phyche_index, all_property, extra_phyche_index)

    return sequence_list, phyche_value




def GetDAC(input_data,**kwargs):
    """
    #################################################################
    Make DAC dictionary.

    :param input_data: file object or sequence list.
    :param phyche_index: physicochemical properties list.
    :param all_property: bool, choose all physicochemical properties or not.
    :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                               It means user-defined phyche_index.
    #################################################################
    """
    
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 2
    
    if 'lag' in kwargs:
        lag = kwargs['lag']
    else:
        lag = 2
    if 'phyche_index' in kwargs:
        phyche_index = kwargs['phyche_index']
    else:
        phyche_index = None
    if 'all_property' in kwargs:
        all_property = kwargs['all_property']
    else:
        all_property = False
    if 'extra_phyche_index' in kwargs:
        extra_phyche_index = kwargs['extra_phyche_index']
    else:
        extra_phyche_index = None
        
    input_data = [input_data]
    sequence_list, phyche_value = ReadyAcc(input_data, k, phyche_index, all_property, extra_phyche_index)
    from PyDNAacutil import MakeACVector
    vec = MakeACVector(sequence_list, lag, phyche_value, k)
    dict_keys = ['DAC_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res



def GetDCC(input_data,**kwargs):
    """
    #################################################################
    Make DCC vector.

    :param input_data: file object or sequence list.
    :param phyche_index: physicochemical properties list.
    :param all_property: bool, choose all physicochemical properties or not.
    :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                               It means user-defined phyche_index.
    #################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 2
    
    if 'lag' in kwargs:
        lag = kwargs['lag']
    else:
        lag = 2
    if 'phyche_index' in kwargs:
        phyche_index = kwargs['phyche_index']
    else:
        phyche_index = None
    if 'all_property' in kwargs:
        all_property = kwargs['all_property']
    else:
        all_property = False
    if 'extra_phyche_index' in kwargs:
        extra_phyche_index = kwargs['extra_phyche_index']
    else:
        extra_phyche_index = None
        
    input_data = [input_data]
    sequence_list, phyche_value = ReadyAcc(input_data, k, phyche_index, all_property, extra_phyche_index)
    from PyDNAacutil import MakeCCVector

    vec = MakeCCVector(sequence_list, lag, phyche_value, k)
    dict_keys = ['DCC_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res
        




def GetDACC(input_data, **kwargs):
    """
    #################################################################
    Make DACC dictionary.

    :param input_data: file object or sequence list.
    :param phyche_index: physicochemical properties list.
    :param all_property: bool, choose all physicochemical properties or not.
    :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                               It means user-defined phyche_index.
    #################################################################
    """
    
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 2
    
    if 'lag' in kwargs:
        lag = kwargs['lag']
    else:
        lag = 2
    if 'phyche_index' in kwargs:
        phyche_index = kwargs['phyche_index']
    else:
        phyche_index = None
    if 'all_property' in kwargs:
        all_property = kwargs['all_property']
    else:
        all_property = False
    if 'extra_phyche_index' in kwargs:
        extra_phyche_index = kwargs['extra_phyche_index']
    else:
        extra_phyche_index = None
        
    input_data = [input_data]
    sequence_list, phyche_value = ReadyAcc(input_data, k, phyche_index, all_property, extra_phyche_index)
    from PyDNAacutil import MakeACVector, MakeCCVector
    zipped = list(zip(MakeACVector(sequence_list, lag, phyche_value, k),
                 MakeCCVector(sequence_list, lag, phyche_value, k)))
    vec = [reduce(lambda x, y: x + y, e) for e in zipped]

    dict_keys = ['DACC_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res



def GetTAC(input_data, **kwargs):
    """
    #################################################################
    Make TAC dictionary.

    :param input_data: file object or sequence list.
    :param phyche_index: physicochemical properties list.
    :param all_property: bool, choose all physicochemical properties or not.
    :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                               It means user-defined phyche_index.
    #################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 3
    
    if 'lag' in kwargs:
        lag = kwargs['lag']
    else:
        lag = 2
    if 'phyche_index' in kwargs:
        phyche_index = kwargs['phyche_index']
    else:
        phyche_index = None
    if 'all_property' in kwargs:
        all_property = kwargs['all_property']
    else:
        all_property = False
    if 'extra_phyche_index' in kwargs:
        extra_phyche_index = kwargs['extra_phyche_index']
    else:
        extra_phyche_index = None
        
    input_data = [input_data]
    sequence_list, phyche_value = ReadyAcc(input_data, k, phyche_index, all_property, extra_phyche_index)
    from PyDNAacutil import MakeACVector
    vec = MakeACVector(sequence_list, lag, phyche_value, k)
    dict_keys = ['TAC_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res



def GetTCC( input_data,**kwargs):
    """
    #################################################################
    Make TCC dictionary.

    :param input_data: file object or sequence list.
    :param phyche_index: physicochemical properties list.
    :param all_property: bool, choose all physicochemical properties or not.
    :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                               It means user-defined phyche_index.
    #################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 3
    
    if 'lag' in kwargs:
        lag = kwargs['lag']
    else:
        lag = 2
    if 'phyche_index' in kwargs:
        phyche_index = kwargs['phyche_index']
    else:
        phyche_index = None
    if 'all_property' in kwargs:
        all_property = kwargs['all_property']
    else:
        all_property = False
    if 'extra_phyche_index' in kwargs:
        extra_phyche_index = kwargs['extra_phyche_index']
    else:
        extra_phyche_index = None
        
    input_data = [input_data]
    sequence_list, phyche_value = ReadyAcc(input_data, k, phyche_index, all_property, extra_phyche_index)
    from PyDNAacutil import MakeCCVector
    vec = MakeCCVector(sequence_list, lag, phyche_value, k)
    dict_keys = ['TCC_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res




def GetTACC(input_data,**kwargs):
    """
    #################################################################
    Make TACC dictionary.

    :param input_data: file object or sequence list.
    :param phyche_index: physicochemical properties list.
    :param all_property: bool, choose all physicochemical properties or not.
    :param extra_phyche_index: dict, the key is the dinucleotide (string), and its corresponding value is a list.
                               It means user-defined phyche_index.
    #################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 3
    
    if 'lag' in kwargs:
        lag = kwargs['lag']
    else:
        lag = 2
    if 'phyche_index' in kwargs:
        phyche_index = kwargs['phyche_index']
    else:
        phyche_index = None
    if 'all_property' in kwargs:
        all_property = kwargs['all_property']
    else:
        all_property = False
    if 'extra_phyche_index' in kwargs:
        extra_phyche_index = kwargs['extra_phyche_index']
    else:
        extra_phyche_index = None
        
    input_data = [input_data]
    sequence_list, phyche_value = ReadyAcc(input_data, k, phyche_index, all_property, extra_phyche_index)

    from PyDNAacutil import MakeACVector, MakeCCVector

    zipped = list(zip(MakeACVector(sequence_list, lag, phyche_value, k),
                 MakeCCVector(sequence_list, lag, phyche_value, k)))
    vec = [reduce(lambda x, y: x + y, e) for e in zipped]
    dict_keys = ['TCC_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res

if __name__ == '__main__':
    extra_phyche_value = {'AA': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11],
                          'AC': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
                          'AG': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
                          'AT': [1.07, 0.22, 0.62, -1.02, 2.51, 1.17],
                          'CA': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
                          'CC': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
                          'CG': [-1.66, -1.22, -0.44, -0.82, -0.29, -1.39],
                          'CT': [0.78, 0.36, 0.09, 0.68, -0.24, -0.62],
                          'GA': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
                          'GC': [-0.08, 0.22, 1.33, -0.35, 0.65, 1.59],
                          'GG': [0.06, 1.08, 0.09, 0.56, -0.82, 0.24],
                          'GT': [1.50, 0.50, 0.80, 0.13, 1.29, 1.04],
                          'TA': [-1.23, -2.37, -0.44, -2.24, -1.51, -1.39],
                          'TC': [-0.08, 0.5, 0.27, 0.13, -0.39, 0.71],
                          'TG': [-1.38, -1.36, -0.27, -0.86, -0.62, -1.25],
                          'TT': [0.06, 0.5, 0.27, 1.59, 0.11, -0.11]}
    phyche_index = \
        [[2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
         [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]

    from PyDNAutil import NormalizeIndex


    dac = GetDAC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'])
    print(dac)
    print(len(dac))

    dac = GetDAC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', all_property=True)
    print(dac)
    print(len(dac))

    dac = GetDAC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'],
                           extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(dac)
    print(len(dac))
    print('\n')


    dcc = GetDCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'])
    print(dcc)
    print(len(dcc))

    dcc = GetDCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',  all_property=True)
    print(dcc)
    print(len(dcc))

    dcc = GetDCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'],
                           extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(dcc)
    print(len(dcc))
    print('\n')

    print('DACC')
    dacc = GetDACC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'])
    print(dacc)
    print(len(dacc))

    dacc = GetDACC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',all_property=True)
    print(dacc)
    print(len(dacc))

    dac = GetDACC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'],
                           extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(dac)
    print(len(dac))
    print('\n')

    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581, 7.176]
    ]

    print('Begin TAC')
    tac = GetTAC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'])
    print(tac)
    print(len(tac))

    tac = GetTAC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',all_property=True)
    print(tac)
    print(len(tac))

    tac = GetTAC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'],
                           extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(tac)
    print(len(tac))
    print('\n')


    print('Begin TCC')
    tcc = GetTCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'])
    print(tcc)
    print(len(tcc))

    tcc = GetTCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',  all_property=True)
    print(tcc)
    print(len(tcc))

    tcc = GetTCC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'],
                           extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(tcc)
    print(len(tcc))
    print('\n')

    print('Begin TACC')
    tacc = GetTACC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'])
    print(tacc)
    print(len(tacc))

    tacc = GetTACC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', all_property=True)
    print(tacc)
    print(len(tacc))

    
    tacc = GetTACC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',  phyche_index=['Dnase I', 'Nucleosome'],
                             extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(tacc)
    print(len(tacc))
    print('\n')