#!/usr/bin/python3.6.7
# -*- coding: utf-8 -*-:v

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	This script read blastn output. It is used one time for the pG4r vs rG4.
        For the other side, see readBlast. The only output is printed on the
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
    bestHit = {} # initialisation of the dictionary
    cptEmpty = 0
    cptHit = 0
    for blastFile in filelist:
        tmp = blastFile.split(':')
        chrQ = tmp[1]
        strandQ = tmp[3].rstrip().split('.')[0]
        namequery = blastFile.split('.')[0]
        locationT = tmp[4]
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
                    chrT = record.alignments[hitN].title.split(':')[0].split('r')[1]
                    locationQ = record.alignments[hitN].title.split(':')[5]
                    if strandT == '+':
                        strandT = '1'
                    else:
                        strandT = '-1'
                    hitN += 1
                    if strandT == strandQ and chrT == chrQ:
                        cptHit += 1
                        hit = True
                        hitLoc = locationQ+'-'+locationT
                        if hitLoc not in dicoHitG4Loc:
                            dicoHitG4Loc[hitLoc] = 0
                        dicoHitG4Loc[hitLoc] += 1
                    # else:
                        # hitLoc = locationQ+'-'+locationT
                        # if hitLoc not in dicoHitLoc:
                            # dicoHitLoc[hitLoc] = 0
                        # dicoHitLoc[hitLoc] += 1
                # if not hit and hitN < maxHit-1:
                    # hitN = 0
                    # while hitN < maxHit:
                        # print(record.alignments[0])
                        # strandT = record.alignments[hitN].title.split(':')[4]
                        # chrT = record.alignments[hitN].title.split(':')[2].split('r')[1]
                        # locationQ = record.alignments[hitN].title.split(':')[6]
                        # if strandT == '+':
                            # strandT = '1'
                        # else:
                            # strandT = '-1'
                        # hitN += 1
            else:
                cptEmpty += 1
    print(nHit)
    print('Hit : ' + str(cptHit))
    print('No hit at all : ' + str(cptEmpty))
    pprint(dicoHitG4Loc)
    #pprint(dicoHitLoc)
    print('--------------')

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    nHit = int(arg.nHit)
    GetBestHit(path, nHit)
