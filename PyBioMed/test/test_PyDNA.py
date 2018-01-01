# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
The script is used for testing.

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.14

Email: gadsby@163.com 
"""
import time
from PyBioMed.PyDNA.PyDNAutil import NormalizeIndex

from PyBioMed.PyDNA.PyDNAac import GetDAC
from PyBioMed.PyDNA.PyDNAac import GetDCC
from PyBioMed.PyDNA.PyDNAac import GetDACC
from PyBioMed.PyDNA.PyDNAac import GetTAC
from PyBioMed.PyDNA.PyDNAac import GetTCC
from PyBioMed.PyDNA.PyDNAac import GetTACC

from PyBioMed.PyDNA.PyDNAnac import GetKmer
from PyBioMed.PyDNA.PyDNAnac import GetRevcKmer
from PyBioMed.PyDNA.PyDNAnac import GetIdKmer

from PyBioMed.PyDNA.PyDNApsenac import GetPseDNC
from PyBioMed.PyDNA.PyDNApsenac import GetPseKNC
from PyBioMed.PyDNA.PyDNApsenac import GetPCPseDNC
from PyBioMed.PyDNA.PyDNApsenac import GetPCPseTNC
from PyBioMed.PyDNA.PyDNApsenac import GetSCPseDNC
from PyBioMed.PyDNA.PyDNApsenac import GetSCPseTNC


def test_pydna():
    
    #os.chdir(os.path.join(os.path.dirname(__file__), os.path.pardir))
    

    #==============================================================================
    # PyDNAac
    #==============================================================================
    
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
    
    
    print '...............................................................'
    print 'testing the GetDAC module'
    
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
    
    print '...............................................................'
    print 'testing the GetDCC module'
    
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
    
    
    print '...............................................................'
    print 'testing the DACC module'
    
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
    
    print '...............................................................'
    print 'testing the GetTAC module'
    
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
    
    
    print '...............................................................'
    print 'testing the GetTCC module'
    
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
    
    print '...............................................................'
    print 'testing the TACC module'
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
    #==============================================================================
    # PyDNAnac
    #==============================================================================
    print '...............................................................'
    print 'testing the GetKmer module'
    
    kmer = GetKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2)
    print(kmer)
     
    kmer = GetKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2,normalize=True)
    print(kmer)
    
    kmer = GetKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2,normalize=True,upto=True)
    print(kmer)
    
    
    print '...............................................................'
    print 'testing the GetRevcKmer module'
    revckmer = GetRevcKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2, normalize=False, upto=False)
    print(revckmer)
    
    
    revckmer = GetRevcKmer('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',k=2,normalize=True, upto=True)
    print(revckmer)
    print('\n')
    
      
    
    
    #==============================================================================
    # PyDNApsenac
    #==============================================================================
    print '...............................................................'
    print 'testing the GetPseDNC module'
    
    psednc = GetPseDNC('ACCCCA',lamada=2, w=0.05)
    print(psednc)
    
    
    print '...............................................................'
    print 'testing the GetPCPseDNC module'
    
    PC_psednc = GetPCPseDNC('ACCCCA', phyche_index=["Tilt", 'Twist', 'Rise', 'Roll', 'Shift', 'Slide'],lamada=2, w=0.05)
    print(PC_psednc)
    
    
    print '...............................................................'
    print 'testing the GetPCPseTNC module'
    
    pc_psetnc = GetPCPseTNC('ACCCCA', phyche_index=['Dnase I', 'Nucleosome'],lamada=2, w=0.05)
    print(pc_psetnc)
    
    print '...............................................................'
    print 'testing the GetSCPseDNC module'
    
    sc_psednc = GetSCPseDNC('ACCCCCA', phyche_index=['Twist', 'Tilt'],lamada=2, w=0.05)
    print(sc_psednc)
    
    print '...............................................................'
    print 'testing the GetSCPseTNC module'
    
    sc_psetnc = GetSCPseTNC('ACCCCCA', phyche_index=['Dnase I', 'Nucleosome'],lamada=1, w=0.05)
    print(sc_psetnc)
    
    sc_psetnc = GetSCPseTNC('ACCCCA', phyche_index=["Dnase I", 'Nucleosome'], lamada=2, w=0.05)
    print(sc_psetnc)
    
    
    phyche_index = [[1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455, -0.070, -0.918, 1.603,
                     -0.654, 0.567, 1.019]]
    
    
    print '...............................................................'
    print 'testing the GetPseDNC module'
    
    dic = GetPseDNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC')
    print(dic)
    print(len(dic))
    
    print '...............................................................'
    print 'testing the GetPseKNC module'
    
    dic = GetPseKNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC')
    print(dic)
    print(len(dic))
    
    print '...............................................................'
    print 'testing the PC-PseDNC module'
    
    print('PC-PseDNC')
    dic = GetPCPseDNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'])
    print(dic)
    print(len(dic))
    
    print '...............................................................'
    print 'testing the GetPCPseTNC module'
    
    dic = GetPCPseTNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC',lamada=1, w=0.05,k=2,phyche_index=['Twist', 'Tilt'])
    print(dic)
    print(len(dic))
    
    
      
    
    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581, 7.176]
    ]
    
    print '...............................................................'
    print 'testing the GetPCPseTNC module'
    
    dic = GetPCPseTNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'],
                                      extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(dic)
    print(len(dic))
    
    
    print '...............................................................'
    print 'testing the GetSCPseDNC module'
    
    dic = GetSCPseDNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'])
    print(dic)
    print(len(dic))
    
    print '...............................................................'
    print 'testing the GetSCPseDNC module'
    
    dic = GetSCPseDNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', all_property=True,lamada=2, w=0.05)
    print(dic)
    print(len(dic))
    
    phyche_index = [[1.019, -0.918, 0.488, 0.567, 0.567, -0.070, -0.579, 0.488, -0.654, -2.455, -0.070, -0.918, 1.603,
                     -0.654, 0.567, 1.019]]
    
    print '...............................................................'
    print 'testing the GetSCPseDNC module'
    
    dic = GetSCPseDNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Twist', 'Tilt'],
                                      extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(dic)
    print(len(dic))
    print()
    
    print '...............................................................'
    print 'testing the GetSCPseTNC module'
    
    dic= GetSCPseTNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'])
    print(dic)
    print(len(dic))
    
    print '...............................................................'
    print 'testing the GetSCPseTNC module'
    
    dic = GetSCPseTNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', all_property=True,lamada=2, w=0.05)
    print(dic)
    print(len(dic))
    
    phyche_index = [
        [7.176, 6.272, 4.736, 7.237, 3.810, 4.156, 4.156, 6.033, 3.410, 3.524, 4.445, 6.033, 1.613, 5.087, 2.169, 7.237,
         3.581, 3.239, 1.668, 2.169, 6.813, 3.868, 5.440, 4.445, 3.810, 4.678, 5.440, 4.156, 2.673, 3.353, 1.668, 4.736,
         4.214, 3.925, 3.353, 5.087, 2.842, 2.448, 4.678, 3.524, 3.581, 2.448, 3.868, 4.156, 3.467, 3.925, 3.239, 6.272,
         2.955, 3.467, 2.673, 1.613, 1.447, 3.581, 3.810, 3.410, 1.447, 2.842, 6.813, 3.810, 2.955, 4.214, 3.581, 7.176]
    ]
    
    dic = GetSCPseTNC('GACTGAACTGCACTTTGGTTTCATATTATTTGCTC', phyche_index=['Dnase I', 'Nucleosome'],
                                      extra_phyche_index=NormalizeIndex(phyche_index, is_convert_dict=True))
    print(dic)
    print(len(dic))
    
    # Normalize PseDNC index Twist, Tilt, Roll, Shift, Slide, Rise.
    original_phyche_value = [
        [0.026, 0.036, 0.031, 0.033, 0.016, 0.026, 0.014, 0.031, 0.025, 0.025, 0.026, 0.036, 0.017, 0.025, 0.016, 0.026],
        [0.038, 0.038, 0.037, 0.036, 0.025, 0.042, 0.026, 0.037, 0.038, 0.036, 0.042, 0.038, 0.018, 0.038, 0.025, 0.038],
        [0.020, 0.023, 0.019, 0.022, 0.017, 0.019, 0.016, 0.019, 0.020, 0.026, 0.019, 0.023, 0.016, 0.020, 0.017, 0.020],
        [1.69, 1.32, 1.46, 1.03, 1.07, 1.43, 1.08, 1.46, 1.32, 1.20, 1.43, 1.32, 0.72, 1.32, 1.07, 1.69],
        [2.26, 3.03, 2.03, 3.83, 1.78, 1.65, 2.00, 2.03, 1.93, 2.61, 1.65, 3.03, 1.20, 1.93, 1.78, 2.26],
        [7.65, 8.93, 7.08, 9.07, 6.38, 8.04, 6.23, 7.08, 8.56, 9.53, 8.04, 8.93, 6.23, 8.56, 6.38, 7.65]]
    for e in NormalizeIndex(original_phyche_value, is_convert_dict=True).items():
        print(e)


if __name__ == '__main__':
    test_pydna()


















