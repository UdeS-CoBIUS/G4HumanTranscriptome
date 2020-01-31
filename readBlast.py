#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	This script read blastn output. It is used one time for the pG4r vs rG4
        blast. For the other side, see readBlastpG4Q. The only output is printed on the
        terminal, and it is : the number of 'perfect hit', and the number of
        blast file without any hit. A perfect hit is when chromosome and strand
        from both blasted fasta are the same.

Usage:
    python ~/PATH/readBlast.py
"""

import sys
import os
import re
import argparse
from pprint import pprint
from Bio.Blast import NCBIXML

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'readBlast')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-nHit', '--nHit', default = 5)
    return parser

def GetBestHit(path, nHit):
    dicoHitG4Loc = {}
    dicoHitLoc= {}
    filelist = sorted(os.listdir(path))# recuperation of all file name in the directory
    cptEmpty = 0
    cptHit = 0
    for blastFile in filelist:
        tmp = blastFile.split(':')
        chrQ = tmp[0].split('r')[1]
        strandQ = tmp[3].rstrip()
        locationQ = blastFile.split(':')[5]
        if strandQ == '+':
            strandQ = '1'
        else:
            strandQ = '-1'
        namequery = blastFile.split('.')[0]
        blastFile = path+blastFile
        blast = NCBIXML.parse(open(blastFile,'rU'))
        for record in blast:
            if record.alignments:
                hit = False
                hitN = 0
                if len(record.alignments) < nHit:
                    maxHit = len(record.alignments) -1
                else:
                    maxHit = nHit
                while not hit and hitN < maxHit:
                    # print(record.alignments[0])
                    strandT = record.alignments[hitN].title.split(':')[3]
                    chrT = record.alignments[hitN].title.split(':')[1]
                    locationT = record.alignments[hitN].title.split(':')[4]
                    hitN += 1
                    if strandT == strandQ and chrT == chrQ:
                        cptHit += 1
                        hit = True
                        hitLoc = locationQ+'-'+locationT
                        if hitLoc not in dicoHitG4Loc:
                            dicoHitG4Loc[hitLoc] = 0
                        dicoHitG4Loc[hitLoc] += 1
            else:
                cptEmpty += 1
    print('Hit : ' + str(cptHit))
    print('No hit at all : ' + str(cptEmpty))
    # pprint(dicoHitG4Loc)

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    nHit = int(arg.nHit)
    GetBestHit(path, nHit)
