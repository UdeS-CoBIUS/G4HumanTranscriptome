#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	August 2018

Description:
	This script aims to predict pG4r in piwi RNA and to computes
		statistiques about it. Firstly, windows from G4RNA screener are filtered
		into pG4r. Then density is computed. The output is printed on the
		screen. Just for information, all piRNA are smaller than 60nt, so there
		is no overlaping windows or pG4r.

Data availability:
	* fasta sequences used with G4RNA screener were downloaded from piRBASE in August 2018.
	* G4RNA screener was launch in August 2018. Results may change if G4RNA is update later.

"""

import sys
import os
import re
import argparse
import recurentFunction as rF

def ReturnG4InGene(inputfile, dicoParam):
    """Merges all windows from windows upper thresholds, it will pG4r.

    This function browses all windows returned by G4RNA Screener and
    keep only those over the thresholds. Here there are no overlaping windows
    because they don't reach 60nt.

    :param inputfile: name of the outputfile of G4RNA screener.
    :type inputfile: string
    :param dicoParam: all parameters given to G4RNA screener.
    :type dicoParam: dictionary

    :returns: G4DetectedPiwi, contains all predicted G4 in genes,
        with there scores.
    :rtype: dictionary
    """
    G4DetectedPiwi = {}
    rnaType = 'piRNA'
    localisation = 'ExonNC'
    oldPassed = False
    passed = False
    inputfile = open(inputfile,"r")
    for line in inputfile:
        if (re.search('^[0-9]', line)):
            words = line.split('\t')
            numRow = words[0]
            infopiRNA = words[1]
            cGcC = float(words[2])
            g4H = float(words[3])
            sequence = words[4]
            startWindow = int(words[5])
            endWindow = int(words[6])
            g4NN = float(words[7].rstrip())
            if (cGcC >= dicoParam['cGcC'] and
                g4H >= dicoParam['g4H'] and
                g4NN >= dicoParam['g4NN']):
                # if window contain value > several threshold
                G4DetectedPiwi[infopiRNA] = [str(cGcC), str(g4H), sequence,
                                            str(g4NN), localisation, rnaType]
    inputfile.close()
    return G4DetectedPiwi

def getLengthpiRNA (inputfile):
    """Counts the length of all piRNA present in the input file.

    :param inputfile: name of the file containing the liste of informations
        for each piRNA.
	:type inputfile: string

    :returns: length, all piRNA length present in the file.
    :rtype: integer
    """
    length = 0
    inputfile = open(inputfile,"r")
    for line in inputfile:
        if (re.search('^[0-9]', line)):
            words = line.split('\t')
            sequence = words[4]
            length = length + len(sequence)
    inputfile.close()
    return length

def getNbrpiRNA(inputfile):
    """Gets the number of all piRNA present in the input file.

    :param inputfile: name of the file containing the list of informations
        for each piRNA.
	:type inputfile: string

    :returns: listepiRNA, list of id piRNA.
    :rtype: list
    """
    listepiRNA = []
    inputfile = open(inputfile,"r")
    for line in inputfile:
        if (re.search('^[0-9]', line)):
            words=line.split('\t')
            infopiRNA = words[1].rstrip()
            listepiRNA.append(infopiRNA)
    inputfile.close()
    return listepiRNA

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'G4AnnotationPIRNA')
    parser.add_argument ('-i', '--inputfile', default = \
        '/home/anais/Documents/Data/Human/G4Conserve-master/' + \
        'piRNA/piR_human_v1.0.csv')
    parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
    parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
    parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
    parser.add_argument ('-E', '--EXTENSION', default = 100)
    parser.add_argument ('-W', '--WINDOW', default = 60)
    parser.add_argument ('-S', '--STEP', default = 10)
    return parser

def main(inputfile, dicoParam):
    piRNA = []
    length = 0
    G4DetectedPiwi = ReturnG4InGene(inputfile, dicoParam)
    length = getLengthpiRNA (inputfile)
    piRNA = getNbrpiRNA(inputfile)
    print '\nInfo PG4 detected in piRNA:'
    print '\tName File  :'+inputfile
    if not G4DetectedPiwi:
        print '        NO PG4 DETECTED'
    print '\nNumber PG4 detected in piRNA:    ' + str(len(G4DetectedPiwi))
    print 'Number piRNA with PG4:    ' + str(len(G4DetectedPiwi))
    print 'Number piRNA tested:    ' + str(len(piRNA))
    print 'Length of all piRNA (nt):    ' + str(length)
    print 'Densité in Kb :    ' \
        + str(round(len(G4DetectedPiwi)/ float(length) * 1000,2))
    print 'Percent piRNA with Kb:    ' + \
        str(round(len(G4DetectedPiwi)/ float(len(piRNA)) * 100,2))

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    inputfile = arg.inputfile
    dicoParam = rF.createDicoParam(arg)
    main(inputfile, dicoParam)
