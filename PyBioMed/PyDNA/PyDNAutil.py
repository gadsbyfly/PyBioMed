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

Date: 2016.07.04

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""


import sys

ALPHABET = 'ACGT'


"""Used for process original data."""


class Seq:
    def __init__(self, name, seq, no):
        self.name = name
        self.seq = seq.upper()
        self.no = no
        self.length = len(seq)

    def __str__(self):
        """Output seq when 'print' method is called."""
        return "%s\tNo:%s\tlength:%s\n%s" % (self.name, str(self.no), str(self.length), self.seq)


def IsUnderAlphabet(s, alphabet):
    """
    #################################################################
    Judge the string is within the scope of the alphabet or not.

    :param s: The string.
    :param alphabet: alphabet.

    Return True or the error character.
    #################################################################
    """
    for e in s:
        if e not in alphabet:
            return e

    return True


def IsFasta(seq):
    """
    #################################################################
    Judge the Seq object is in FASTA format.
    Two situation:
    1. No seq name.
    2. Seq name is illegal.
    3. No sequence.

    :param seq: Seq object.
    #################################################################
    """
    if not seq.name:
        error_info = 'Error, sequence ' + str(seq.no) + ' has no sequence name.'
        print(seq)
        sys.stderr.write(error_info)
        return False
    if -1 != seq.name.find('>'):
        error_info = 'Error, sequence ' + str(seq.no) + ' name has > character.'
        sys.stderr.write(error_info)
        return False
    if 0 == seq.length:
        error_info = 'Error, sequence ' + str(seq.no) + ' is null.'
        sys.stderr.write(error_info)
        return False

    return True


def ReadFasta(f):
    """
    #################################################################
    Read a fasta file.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return Seq obj list.
    #################################################################
    """
    name, seq = '', ''
    count = 0
    seq_list = []
    lines = f.readlines()
    for line in lines:
        if not line:
            break

        if '>' == line[0]:
            if 0 != count or (0 == count and seq != ''):
                if IsFasta(Seq(name, seq, count)):
                    seq_list.append(Seq(name, seq, count))
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    count += 1
    if IsFasta(Seq(name, seq, count)):
        seq_list.append(Seq(name, seq, count))
    else:
        sys.exit(0)

    return seq_list


def ReadFastaYield(f):
    """
    #################################################################
    Yields a Seq object.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)
    #################################################################
    """
    name, seq = '', ''
    count = 0
    while True:
        line = f.readline()
        if not line:
            break

        if '>' == line[0]:
            if 0 != count or (0 == count and seq != ''):
                if IsFasta(Seq(name, seq, count)):
                    yield Seq(name, seq, count)
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    if IsFasta(Seq(name, seq, count)):
        yield Seq(name, seq, count)
    else:
        sys.exit(0)


def ReadFastaCheckDna(f):
    """
    #################################################################
    Read the fasta file, and check its legality.

    :param f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return the seq list.
    #################################################################
    """
    seq_list = []
    for e in ReadFastaYield(f):
        # print e
        res = IsUnderAlphabet(e.seq, ALPHABET)
        if res:
            seq_list.append(e)
        else:
            error_info = 'Sorry, sequence ' + str(e.no) \
                         + ' has character ' + str(res) + '.(The character must be A or C or G or T)'
            sys.stderr(error_info)
            sys.exit(0)

    return seq_list


def GetSequenceCheckDna(f):
    """
    #################################################################
    Read the fasta file.

    Input: f: HANDLE to input. e.g. sys.stdin, or open(<file>)

    Return the sequence list.
    #################################################################
    """
    sequence_list = []
    for e in ReadFastaYield(f):
        # print e
        res = IsUnderAlphabet(e.seq, ALPHABET)
        if res is not True:
            error_info = 'Sorry, sequence ' + str(e.no) \
                         + ' has character ' + str(res) + '.(The character must be A, C, G or T)'
            sys.stderr.write(error_info)
            sys.exit(0)
        else:
            sequence_list.append(e.seq)

    return sequence_list


def IsSequenceList(sequence_list):
    """
    #################################################################
    Judge the sequence list is within the scope of alphabet and 
    change the lowercase to capital.
    #################################################################
    """
    count = 0
    new_sequence_list = []

    for e in sequence_list:
        e = e.upper()
        count += 1
        res = IsUnderAlphabet(e, ALPHABET)
        if res is not True:
            error_info = 'Sorry, sequence ' + str(count) \
                         + ' has illegal character ' + str(res) + '.(The character must be A, C, G or T)'
            sys.stderr.write(error_info)
            return False
        else:
            new_sequence_list.append(e)

    return new_sequence_list


def GetData(input_data, desc=False):
    """
    #################################################################
    Get sequence data from file or list with check.

    :param input_data: type file or list
    :param desc: with this option, the return value will be a Seq object list(it only works in file object).
    :return: sequence data or shutdown.
    #################################################################
    """
    if hasattr(input_data, 'read'):
        if desc is False:
            return GetSequenceCheckDna(input_data)
        else:
            return ReadFastaCheckDna(input_data)
    elif isinstance(input_data, list):
        input_data = IsSequenceList(input_data)
        if input_data is not False:
            return input_data
        else:
            sys.exit(0)
    else:
        error_info = 'Sorry, the parameter in get_data method must be list or file type.'
        sys.stderr.write(error_info)
        sys.exit(0)


"""Some basic function for generate feature vector."""


def Frequency(tol_str, tar_str):
    """
    #################################################################
    Generate the frequency of tar_str in tol_str.

    :param tol_str: mother string.
    :param tar_str: substring.
    #################################################################
    """
    i, j, tar_count = 0, 0, 0
    len_tol_str = len(tol_str)
    len_tar_str = len(tar_str)
    while i < len_tol_str and j < len_tar_str:
        if tol_str[i] == tar_str[j]:
            i += 1
            j += 1
            if j >= len_tar_str:
                tar_count += 1
                i = i - j + 1
                j = 0
        else:
            i = i - j + 1
            j = 0

    return tar_count


def WriteLibsvm(vector_list, label_list, write_file):
    """
    #################################################################
    Write the vector into disk in libSVM format.
    #################################################################
    """
    len_vector_list = len(vector_list)
    len_label_list = len(label_list)
    if len_vector_list == 0:
        sys.stderr.write("The vector is none.")
        sys.exit(1)
    if len_label_list == 0:
        sys.stderr.write("The label is none.")
        sys.exit(1)
    if len_vector_list != len_label_list:
        sys.stderr.write("The length of vector and label is different.")
        sys.exit(1)

    with open(write_file, 'w') as f:
        len_vector = len(vector_list[0])
        for i in range(len_vector_list):
            temp_write = str(label_list[i])
            for j in range(0, len_vector):
                temp_write += ' ' + str(j + 1) + ':' + str(vector_list[i][j])
            f.write(temp_write)
            f.write('\n')


def GeneratePhycheValue(k, phyche_index=None, all_property=False, extra_phyche_index=None):
    """
    #################################################################
    Combine the user selected phyche_list, is_all_property and 
    extra_phyche_index to a new standard phyche_value.
    #################################################################
    """
    if phyche_index is None:
        phyche_index = []
    if extra_phyche_index is None:
        extra_phyche_index = {}

    diphyche_list = ['Base stacking', 'Protein induced deformability', 'B-DNA twist', 'Dinucleotide GC Content',
                     'A-philicity', 'Propeller twist', 'Duplex stability:(freeenergy)',
                     'Duplex tability(disruptenergy)', 'DNA denaturation', 'Bending stiffness', 'Protein DNA twist',
                     'Stabilising energy of Z-DNA', 'Aida_BA_transition', 'Breslauer_dG', 'Breslauer_dH',
                     'Breslauer_dS', 'Electron_interaction', 'Hartman_trans_free_energy', 'Helix-Coil_transition',
                     'Ivanov_BA_transition', 'Lisser_BZ_transition', 'Polar_interaction', 'SantaLucia_dG',
                     'SantaLucia_dH', 'SantaLucia_dS', 'Sarai_flexibility', 'Stability', 'Stacking_energy',
                     'Sugimoto_dG', 'Sugimoto_dH', 'Sugimoto_dS', 'Watson-Crick_interaction', 'Twist', 'Tilt',
                     'Roll', 'Shift', 'Slide', 'Rise']
    triphyche_list = ['Dnase I', 'Bendability (DNAse)', 'Bendability (consensus)', 'Trinucleotide GC Content',
                      'Nucleosome positioning', 'Consensus_roll', 'Consensus-Rigid', 'Dnase I-Rigid', 'MW-Daltons',
                      'MW-kg', 'Nucleosome', 'Nucleosome-Rigid']

    # Set and check physicochemical properties.
    if 2 == k:
        if all_property is True:
            phyche_index = diphyche_list
        else:
            for e in phyche_index:
                if e not in diphyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    import sys

                    sys.stderr.write(error_info)
                    sys.exit(0)
    elif 3 == k:
        if all_property is True:
            phyche_index = triphyche_list
        else:
            for e in phyche_index:
                if e not in triphyche_list:
                    error_info = 'Sorry, the physicochemical properties ' + e + ' is not exit.'
                    import sys

                    sys.stderr.write(error_info)
                    sys.exit(0)

    # Generate phyche_value.
    from PyDNApsenacutil import GetPhycheIndex, ExtendPhycheIndex

    return ExtendPhycheIndex(GetPhycheIndex(k, phyche_index), extra_phyche_index)


def ConvertPhycheIndexToDict(phyche_index):
    """
    #################################################################
    Convert phyche index from list to dict.
    #################################################################
    """
    # for e in phyche_index:
    #     print e
    len_index_value = len(phyche_index[0])
    k = 0
    for i in range(1, 10):
        if len_index_value < 4**i:
            error_infor = 'Sorry, the number of each index value is must be 4^k.'
            sys.stdout.write(error_infor)
            sys.exit(0)
        if len_index_value == 4**i:
            k = i
            break
    from PyDNAnacutil import MakeKmerList
    kmer_list = MakeKmerList(k, ALPHABET)
    # print kmer_list
    len_kmer = len(kmer_list)
    phyche_index_dict = {}
    for kmer in kmer_list:
        phyche_index_dict[kmer] = []
    # print phyche_index_dict
    phyche_index = list(zip(*phyche_index))
    for i in range(len_kmer):
        phyche_index_dict[kmer_list[i]] = list(phyche_index[i])

    return phyche_index_dict


def StandardDeviation(value_list):
    """
    #################################################################
    Return standard deviation.
    #################################################################
    """
    from math import sqrt
    from math import pow
    n = len(value_list)
    average_value = sum(value_list) * 1.0 / n
    return sqrt(sum([pow(e - average_value, 2) for e in value_list]) * 1.0 / (n-1))


def NormalizeIndex(phyche_index, is_convert_dict=False):
    """
    #################################################################
    Normalize the physicochemical index.
    #################################################################
    """
    normalize_phyche_value = []
    for phyche_value in phyche_index:
        average_phyche_value = sum(phyche_value) * 1.0 / len(phyche_value)
        sd_phyche = StandardDeviation(phyche_value)
        normalize_phyche_value.append([round((e - average_phyche_value) / sd_phyche, 2) for e in phyche_value])

    if is_convert_dict is True:
        return ConvertPhycheIndexToDict(normalize_phyche_value)

    return normalize_phyche_value

if __name__ == "__main__":

    re = ['GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', 'GACTGAACTGCACTTTGGTTTCATATTATTTGCTC']


    phyche_index = [[1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455,-0.070, -0.918, 1.603, -0.654, 0.567, 1.019]]
    print (NormalizeIndex(phyche_index,is_convert_dict = False)[0])
    import os    
    print os.path.abspath('.')
    
