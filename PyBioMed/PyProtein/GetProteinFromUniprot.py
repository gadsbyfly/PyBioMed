# -*- coding: utf-8 -*-
"""
################################################################################################

This module is used to download the protein sequence from the uniprot (http://www.uniprot.org/) 

website. You can only need input a protein ID or prepare a file (ID.txt) related to ID. You can

 obtain a .txt (ProteinSequence.txt) file saving protein sequence you need.  You can freely use 
 
 and distribute it. If you hava  any problem, you could contact with us timely!
 
Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com

################################################################################################
"""

import urllib
import string


##################################################################################################
def GetProteinSequence(ProteinID):
    """
    #########################################################################################
    Get the protein sequence from the uniprot website by ID.

    Usage:

    result=GetProteinSequence(ProteinID)

    Input: ProteinID is a string indicating ID such as "P48039".

    Output: result is a protein sequence.
    #########################################################################################
    """

    ID = str(ProteinID)
    localfile = urllib.urlopen('http://www.uniprot.org/uniprot/' + ID + '.fasta')
    temp = localfile.readlines()
    res = ''
    for i in range(1, len(temp)):
        res = res + string.strip(temp[i])
    return res


##################################################################################################
def GetProteinSequenceFromTxt(path, openfile, savefile):
    """
    #########################################################################################
    Get the protein sequence from the uniprot website by the file containing ID.

    Usage:

    result=GetProteinSequenceFromTxt(path,openfile,savefile)

    Input: path is a directory path containing the ID file such as "/home/orient/protein/"

    openfile is the ID file such as "proteinID.txt"

    savefile is the file saving the obtained protein sequences such as "protein.txt"
    #########################################################################################
    """
    f1 = file(path + savefile, 'wb')
    f2 = file(path + openfile, 'r')
    #	res=[]
    for index, i in enumerate(f2):

        itrim = string.strip(i)
        if itrim == "":
            continue
        else:
            temp = GetProteinSequence(itrim)
            print "--------------------------------------------------------"
            print "The %d protein sequence has been downloaded!" % (index + 1)
            print temp
            f1.write(temp + '\n')
            print "--------------------------------------------------------"
        #		res.append(temp+'\n')
        #	f1.writelines(res)
    f2.close()
    f1.close()
    return 0


##################################################################################################
if __name__ == '__main__':

    localfile = ['P48039']
    for index, i in enumerate(localfile):
        itrim = string.strip(i)
        if itrim == "":
            continue
        else:
            temp = GetProteinSequence(itrim)
            print "--------------------------------------------------------"
            print "The %d protein sequence has been downloaded!" % (index + 1)
            print temp
            print "--------------------------------------------------------"
