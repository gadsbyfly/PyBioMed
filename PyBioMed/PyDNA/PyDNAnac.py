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

Date: 2016.10.11

Email: gadsby@163.com and oriental-cds@163.com

##############################################################################
"""



from PyDNAnacutil import MakeUptoKmerList, MakeRevcompKmerList, MakeKmerVector
from PyDNAutil import GetData


def CheckNacPara(k, normalize=False, upto=False, alphabet='ACGT'):
    """
    ###########################################################################
    Check the nac parameter's validation.
    ###########################################################################
    """
    try:
        if not isinstance(k, int) or k <= 0:
            raise ValueError("Error, parameter k must be an integer and larger than 0.")
        elif not isinstance(normalize, bool):
            raise ValueError("Error, parameter normalize must be bool type.")
        elif not isinstance(upto, bool):
            raise ValueError("Error, parameter upto must be bool type.")
        elif alphabet != 'ACGT':
            raise ValueError("Error, parameter alphabet must be 'ACGT'.")
    except ValueError:
        raise


def GetKmerList(k, upto, alphabet):
    """
    ###########################################################################
    Get the kmer list.

    :param k: int, the k value of kmer, it should be larger than 0.
    :param upto: bool, whether to generate 1-kmer, 2-kmer, ..., k-mer.
    :param alphabet: string.
    ###########################################################################
    """
    if upto:
        k_list = list(range(1, k + 1))
    else:
        k_list = list(range(k, k + 1))
    kmer_list = MakeUptoKmerList(k_list, alphabet)

    return kmer_list


def GetKmer(data, **kwargs):
    """
    ###########################################################################
    Make a kmer dictionary with options k, upto, revcomp, normalize.
    
    :param k: int, the k value of kmer, it should be larger than 0.
    :param normalize: bool, normalize the result vector or not.
    :param upto: bool, whether to generate 1-kmer, 2-kmer, ..., k-mer.
    :param alphabet: string.
    :param data: file object or sequence list.
    :return: kmer vector.
    ###########################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 1            
    if 'normalize' in kwargs:
        normalize = kwargs['normalize']
    else:
        normalize = False
    if 'upto' in kwargs:
        upto =kwargs['upto']
    else:
        upto = False        
    if 'alphabet' in kwargs:
        alphabet = kwargs['alphabet']
    else:
        alphabet = "ACGT"
    
    data = [data]    
    sequence_list = GetData(data)

    kmer_list = GetKmerList(k, upto, alphabet)

    rev_kmer_list = []
    revcomp = False
    vec = MakeKmerVector(sequence_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
    
    dict_keys = ['Kmer_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res


def GetRevcKmer(data, **kwargs):
    """
    ###########################################################################
    Make a reverse compliment kmer dictionary with options k, upto, normalize.

    :param data: file object or sequence list.
    :return: reverse compliment kmer vector.
    ###########################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 1            
    if 'normalize' in kwargs:
        normalize = kwargs['normalize']
    else:
        normalize = False
    if 'upto' in kwargs:
        upto =kwargs['upto']
    else:
        upto = False        
    if 'alphabet' in kwargs:
        alphabet = kwargs['alphabet']
    else:
        alphabet = "ACGT"
        
    
    data = [data]
    sequence_list = GetData(data)

    kmer_list = GetKmerList(k, upto, alphabet)

    # Use lexicographically first version of {kmer, revcomp(kmer)}.
    rev_kmer_list = MakeRevcompKmerList(kmer_list)
    revcomp = True
    vec = MakeKmerVector(sequence_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
    
    dict_keys = ['RevcKmer_%s'%i for i in range(1,len(vec[0])+1)]
    res = dict(zip(dict_keys,vec[0]))
    return res




def GetIdKmer(data, hs, non_hs,**kwargs):
    """
    ###########################################################################
    Make IDKmer vector.

    :param data: Need to processed FASTA file.
    :param hs: Positive FASTA file.
    :param non_hs: Negative FASTA file.
    :param k: int, the k value of kmer, it should be larger than 0.
    :param upto: bool, whether to generate 1-kmer, 2-kmer, ..., k-mer.
    :param alphabet: string.
    ###########################################################################
    """
    if 'k' in kwargs:
        k = kwargs['k']
    else:
        k = 6            
    if 'upto' in kwargs:
        upto =kwargs['upto']
    else:
        upto = True        
    if 'alphabet' in kwargs:
        alphabet = kwargs['alphabet']
    else:
        alphabet = "ACGT"
        
    
    from PyDNAnacutil import MakeKmerList
    from PyDNAnacutil import Diversity
    from PyDNAnacutil import IdXS

    rev_kmer_list, upto, revcomp, normalize = [], False, False, False

    pos_s_list = GetData(hs)
    neg_s_list = GetData(non_hs)
    # print k
    if upto is False:
        k_list = [k]
    else:
        k_list = list(range(1, k+1))

    # print 'k_list =', k_list

    # Get all kmer ID from 1-kmer to 6-kmer.
    # Calculate standard source S vector.
    pos_s_vec, neg_s_vec = [], []
    diversity_pos_s, diversity_neg_s = [], []
    for k in k_list:
        kmer_list = MakeKmerList(k, alphabet)

        temp_pos_s_vec = MakeKmerVector(pos_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
        temp_neg_s_vec = MakeKmerVector(neg_s_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)

        temp_pos_s_vec = [sum(e) for e in zip(*[e for e in temp_pos_s_vec])]
        temp_neg_s_vec = [sum(e) for e in zip(*[e for e in temp_neg_s_vec])]

        pos_s_vec.append(temp_pos_s_vec)
        neg_s_vec.append(temp_neg_s_vec)

        diversity_pos_s.append(Diversity(temp_pos_s_vec))
        diversity_neg_s.append(Diversity(temp_neg_s_vec))

    # Calculate Diversity(X) and ID(X, S).
    sequence_list = GetData(data)
    vec = []

    for seq in sequence_list:
        # print seq
        temp_vec = []
        for k in k_list:
            kmer_list = MakeKmerList(k, alphabet)
            seq_list = [seq]
            kmer_vec = MakeKmerVector(seq_list, kmer_list, rev_kmer_list, k, upto, revcomp, normalize)
            # print 'k', k
            # print 'kmer_vec', kmer_vec

            # print diversity_pos_s
            if upto is False:
                k = 1

            # print 'pos_vec', pos_s_vec
            # print 'neg_vec', neg_s_vec
            # print 'diversity_pos_s', diversity_pos_s

            temp_vec.append(round(IdXS(kmer_vec[0], pos_s_vec[k-1], diversity_pos_s[k-1]), 3))
            temp_vec.append(round(IdXS(kmer_vec[0], neg_s_vec[k-1], diversity_neg_s[k-1]), 3))

        vec.append(temp_vec)

    return vec


if __name__ == '__main__':
    # kmer =Kmer(k=1)
    # kmer =RevcKmer(k=1, normalize=True, alphabet='ACGT')
    # kmer =IDkmer(k=1)
    
    kmer = GetKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2)
    print(kmer)
 
    kmer = GetKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2,normalize=True)
    print(kmer)

    kmer = GetKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2,normalize=True,upto=True)
    print(kmer)


    revckmer = GetRevcKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2, normalize=False, upto=False)
    print(revckmer)


    revckmer = GetRevcKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2,normalize=True, upto=True)
    print(revckmer)
    print('\n')

