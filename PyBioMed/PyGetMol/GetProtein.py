# -*- coding: utf-8 -*-
#  Copyright (c) 2016-2017, Zhijiang Yao, Jie Dong and Dongsheng Cao
#  All rights reserved.
#  This file is part of the PyBioMed.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the PyBioMed source tree.
"""
Created on Sat Jul 13 11:18:26 2013

This module is used for downloading the PDB file from RCSB PDB web and 

extract its amino acid sequence. This module can also download the protein 

sequence from the uniprot (http://www.uniprot.org/) website. You can only 

need input a protein ID or prepare a file (ID.txt) related to ID. You can

obtain a .txt (ProteinSequence.txt) file saving protein sequence you need.  

Authors: Zhijiang Yao and Dongsheng Cao.

Date: 2016.06.04

Email: gadsby@163.com
"""

import os, ftplib, gzip, sys
import time
from threading import Thread

ThreadStop = Thread._Thread__stop
class TimeoutException(Exception):
    pass
HOSTNAME = "ftp.wwpdb.org"
DIRECTORY = "/pub/pdb/data/structures/all/pdb/"
PREFIX = "pdb"
SUFFIX = ".ent.gz"

# List of three and one letter amino acid codes
_aa_index = [('ALA', 'A'),
             ('CYS', 'C'),
             ('ASP', 'D'),
             ('GLU', 'E'),
             ('PHE', 'F'),
             ('GLY', 'G'),
             ('HIS', 'H'),
             ('HSE', 'H'),
             ('HSD', 'H'),
             ('ILE', 'I'),
             ('LYS', 'K'),
             ('LEU', 'L'),
             ('MET', 'M'),
             ('MSE', 'M'),
             ('ASN', 'N'),
             ('PRO', 'P'),
             ('GLN', 'Q'),
             ('ARG', 'R'),
             ('SER', 'S'),
             ('THR', 'T'),
             ('VAL', 'V'),
             ('TRP', 'W'),
             ('TYR', 'Y')]


# AA3_TO_AA1 = dict(_aa_index)
# AA1_TO_AA3 = dict([(aa[1],aa[0]) for aa in _aa_index])


def timelimited(timeout):
    """
    A decorator to limit the execution time.
    """
    def decorator(function):
        def decorator2(*args,**kwargs):
            class TimeLimited(Thread):
                def __init__(self,_error= None,):
                    Thread.__init__(self)
                    self._error = _error

                def run(self):
                    try:
                        self.result = function(*args,**kwargs)
                    except Exception as e:
                        self._error = e

                def _stop(self):
                    if self.isAlive():
                        ThreadStop(self)

            t = TimeLimited()
            t.start()
            t.join(timeout)

            if isinstance(t._error,TimeoutException):
                t._stop()
                print 'Network timeout, skip!!'
                return 'Network timeout, skip!!'

            if t.isAlive():
                t._stop()
                print 'Network timeout, skip!!'
                return 'Network timeout, skip!!'

            if t._error is None:
                return t.result

        return decorator2
    return decorator


def unZip(some_file, some_output):
    """
    Unzip some_file using the gzip library and write to some_output.
    """

    f = gzip.open(some_file, 'r')
    g = open(some_output, 'w')
    g.writelines(f.readlines())
    f.close()
    g.close()

    os.remove(some_file)


def pdbDownload(file_list, hostname=HOSTNAME, directory=DIRECTORY, prefix=PREFIX,
                suffix=SUFFIX):
    """
    Download all pdb files in file_list and unzip them.
    """

    success = True

    # Log into server
    print "Connecting..."
    ftp = ftplib.FTP()
    ftp.connect(hostname)
    ftp.login()

    # Remove .pdb extensions from file_list
    for file_index, file in enumerate(file_list):
        try:
            file_list[file_index] = file[:file.index(".pdb")]
        except ValueError:
            pass

    # Download all files in file_list
    to_get = ["%s/%s%s%s" % (directory, prefix, f, suffix) for f in file_list]
    to_write = ["%s%s" % (f, suffix) for f in file_list]

    for i in range(len(to_get)):
        try:
            ftp.retrbinary("RETR %s" % to_get[i], open(to_write[i], "wb").write)
            final_name = "%s.pdb" % to_write[i][:to_write[i].index(".")]
            unZip(to_write[i], final_name)
            print "%s retrieved successfully." % final_name
        except ftplib.error_perm:
            os.remove(to_write[i])
            print "ERROR!  %s could not be retrieved!" % file_list[i]
            success = False

    # Log out
    ftp.quit()

    if success:
        return True
    else:
        return False


def GetPDB(pdbidlist=[]):
    """
    Download the PDB file from PDB FTP server by providing a list of pdb id.
    """

    for i in pdbidlist:
        templist = []
        templist.append(i)
        pdbDownload(templist)

    return True


def pdbSeq(pdb, use_atoms=False):
    """
    Parse the SEQRES entries in a pdb file.  If this fails, use the ATOM
    entries.  Return dictionary of sequences keyed to chain and type of
    sequence used.
    """

    # Try using SEQRES
    seq = [l for l in pdb if l[0:6] == "SEQRES"]
    if len(seq) != 0 and not use_atoms:
        seq_type = "SEQRES"
        chain_dict = dict([(l[11], []) for l in seq])
        for c in chain_dict.keys():
            chain_seq = [l[19:70].split() for l in seq if l[11] == c]
            for x in chain_seq:
                chain_dict[c].extend(x)

    # Otherwise, use ATOM
    else:

        seq_type = "ATOM  "

        # Check to see if there are multiple models.  If there are, only look
        # at the first model.
        models = [i for i, l in enumerate(pdb) if l.startswith("MODEL")]
        if len(models) > 1:
            pdb = pdb[models[0]:models[1]]

            # Grab all CA from ATOM entries, as well as MSE from HETATM
        atoms = []
        for l in pdb:
            if l[0:6] == "ATOM  " and l[13:16] == "CA ":

                # Check to see if this is a second conformation of the previous
                # atom
                if len(atoms) != 0:
                    if atoms[-1][17:26] == l[17:26]:
                        continue

                atoms.append(l)
            elif l[0:6] == "HETATM" and l[13:16] == "CA " and l[17:20] == "MSE":

                # Check to see if this is a second conformation of the previous
                # atom
                if len(atoms) != 0:
                    if atoms[-1][17:26] == l[17:26]:
                        continue

                atoms.append(l)

        chain_dict = dict([(l[21], []) for l in atoms])
        for c in chain_dict.keys():
            chain_dict[c] = [l[17:20] for l in atoms if l[21] == c]

    AA3_TO_AA1 = dict(_aa_index)
    tempchain = chain_dict.keys()
    tempchain.sort()
    seq = ''
    for i in tempchain:
        for j in chain_dict[i]:
            if j in AA3_TO_AA1.keys():
                seq = seq + (AA3_TO_AA1[j])
            else:
                seq = seq + ('X')

    return seq, seq_type


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


def GetSeqFromPDB(pdbfile=''):
    """
    Get the amino acids sequences from pdb file.
    """
    f1 = file(pdbfile, 'r')
    res1, res2 = pdbSeq(f1)
    f1.close()
    return res1


class Seq:
    def __init__(self, name, seq, no):
        self.name = name
        self.seq = seq.upper()
        self.no = no
        self.length = len(seq)

    def __str__(self):
        """Output seq when 'print' method is called."""
        return "%s\tNo:%s\tlength:%s\n%s" % (self.name, str(self.no), str(self.length), self.seq)


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
                    seq_list.append(seq)
                else:
                    sys.exit(0)

            seq = ''
            name = line[1:].strip()
            count += 1
        else:
            seq += line.strip()

    count += 1
    if IsFasta(Seq(name, seq, count)):
        seq_list.append(seq)
    else:
        sys.exit(0)

    return seq_list


####################################################################
if __name__ == "__main__":

    print '-'*10+'START'+'-'*10
    @timelimited(10)
    def run_GetPDB():

        GetPDB(['1atp', '1efz', '1f88'])

    @timelimited(10)
    def run_GetSeqFromPDB():

        seq = GetSeqFromPDB('1atp.pdb')
        print seq

    @timelimited(10)
    def run_GetProteinSequence():
        GetProteinSequence('O00560')

        print ReadFasta(open('../test/test_data/protein.fasta'))

    run_GetPDB()
    print '-'*25
    run_GetSeqFromPDB()
    print '-'*25
    run_GetProteinSequence()
    print '-'*10+'END'+'-'*10
