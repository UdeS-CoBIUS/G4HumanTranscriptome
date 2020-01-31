#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	August 2018

Description:
	This script aims to predict pG4r in transfer RNA and to computes
		statistiques about it. Firstly, windows from G4RNA screener are merged
		into pG4r. Then density is computed. The output is printed on the
		screen. Data are dowloaded from Genome tRNA Database in August 2018.

Data availability:
	* fasta sequences used with G4RNA screener were downloaded from Genomic tRNA Database in August 2018.
	* G4RNA screener was launch in August 2018. Results may change if G4RNA is update later.

"""

import sys
import os
import re
import argparse
import recurentFunction as rF

def readLinetRNA(line):
    words = line.split('\t')
    dicoL = {'tRNAInfo' : words[1],
            'cGcC' : float(words[2]),
            'g4H' : float(words[3]),
            'Sequence' : words[4],
            'wStart' : int(words[5]),
            'wEnd' : int(words[6]),
            'g4NN' : float(words[7].rstrip())}
    infoSplit = dicoL['tRNAInfo'].split('(')
    positions = dicoL['tRNAInfo'].split('chr')[2].split(':')
    dicoL.update({'geneId' : infoSplit[0].replace(" ",""),
                'Chromosome' : positions[0],
                'tRNAStart' : int(positions[1].split('-')[0]),
                'tRNAEnd' : int(positions[1].split('(')[0].split('-')[1]),
                'Strand' : dicoL['tRNAInfo'].split('chr')[2].split('(')[1]})
    dicoL['Description'] = dicoL['Chromosome'] + '|' + str(dicoL['tRNAStart'])\
        + '-' + str(dicoL['tRNAEnd'])
    return dicoL

def ReturnG4InGene(inputfile, parametersTool): ## for G4 in gene
    """Merges all windows from windows upper thresholds, it will pG4r.

    This function browsea all windows returned by G4RNA Screener and
    keep only those over the thresholds. Successive windows are merged.

    :param inputfile: name of the outputfile of G4RNA screener.
    :type inputfile: string
    :param dicoParam: all parameters given to G4RNA screener.
    :type dicoParam: dictionary

    :returns: G4DetectedIntRNA, contains all predicted G4 in genes,
        with there scores.
    :rtype: dictionary
    """
    G4DetectedIntRNA = {}
    rnaType = 'tRNA'
    location = 'ExonNC'
    oldPassed = False
    passed = False
    inputfile = open(inputfile,"r")
    for line in inputfile:
        if (re.search('^[0-9]', line)):
            dL = readLinetRNA(line)
            if (dL['cGcC'] >= dicoParam['cGcC'] and
                dL['g4H'] >= dicoParam['g4H'] and
                dL['g4NN'] >= dicoParam['g4NN']):
                passed = True
                if oldPassed != passed:
                    # if it's the first windows
                    descriptionOverThreshold = dL['Description']
                    gene_startG4 = dL['geneId']
                    sequenceG4 = dL['Sequence']
                    oldPassed = passed
                    if dL['Strand'] == '+':
                        startG4 = dL['wStart']
                        endG4 = dL['wEnd']
                    else:
                        startG4 = dL['tRNAEnd'] - \
                            (dL['wStart'] - dL['tRNAStart'])
                        endG4 = startG4 - len(dL['Sequence']) + 1
                    listeCGcC = [ dL['cGcC'] ]
                    listeG4Hunter = [ dL['g4H'] ]
                    listeG4NN = [ dL['g4NN'] ]
                else:
                    # if it's not the first windows above the thresholds
                    sequenceG4 = rF.addWindowToG4Seq(sequenceG4,
                            dL['Sequence'], dicoParam['Step'],
                            dicoParam['Window'])
                    if strand == '+':
                        endG4 = dL['wEnd']
                    else:
                        endG4 = startG4 - len(sequenceG4) + 1
                    listeCGcC.append(dL['cGcC'])
                    listeG4Hunter.append(dL['g4H'])
                    listeG4NN.append(dL['g4NN'])
            if (dL['cGcC'] < dicoParam['cGcC'] or
                dL['g4H'] < dicoParam['g4H'] or
                dL['g4NN'] < dicoParam['g4NN'] or
                descriptionOverThreshold != dL['Description']):
                # one of the score is under his threshold
                # or this windows is from another gene
                passed = False
                if oldPassed != passed:
                    # last windows before under the thresolds
                    meanCGcC = mean(listeCGcC)
                    meanG4Hunter = mean(listeG4Hunter)
                    meanG4NN = mean(listeG4NN)
                    oldPassed = passed
                    headerG4 = rF.createIdG4(dL['geneId'],
                            startG4, endG4, dL['Strand'])
                    if (G4DetectedIntRNA.has_key(headerG4) == False and
                        strand != None):
                        G4DetectedIntRNA[headerG4] = [str(meanCGcC),
                                                    str(meanG4Hunter),
                                                    sequenceG4,
                                                    str(meanG4NN),
                                                    location, rnaType]
    inputfile.close()
    return G4DetectedIntRNA

def getLengthTRNA (inputfile):
    """Counts the length of all tRNA present in the input file.

    :param inputfile: name of the file containing the list of informations
        for each tRNA.
    :type inputfile: string

    :returns: length, length of all tRNA present in the file.
    :rtype: integer
    """
    length = 0
    inputfile = open(inputfile,"r")
    for line in inputfile:
        if (re.search('^[0-9]', line)):
            words = line.split('\t')
            infoTRNA = words[1]
            tRNAStart = int(infoTRNA.split('chr')[2] \
                .split(':')[1].split('-')[0])
            tRNAEnd = int(infoTRNA.split('chr')[2] \
                .split(':')[1].split('(')[0].split('-')[1])
            length += rF.GetLengthFraction(tRNAStart, tRNAEnd)
    inputfile.close()
    return length

def getNbrTRNAWithPG4(dictionary):
    """Counts the number tRNA containing pG4r.

    :param dictionary: informations of region G4 detected for each
        G4 region predicted in tRNA.
    :type dictionary: dictionary

    :returns: nbrTRNAWithPG4, number of tRNA containing pG4r.
    :rtype: list
    """
    listeTRNA = []
    listePG4 = list(dictionary)
    for headerPG4 in listePG4:
        infoTRNA = headerPG4.split('|')[0]
        if infoTRNA not in listeTRNA:
            listeTRNA.append(infoTRNA)
    return listeTRNA

def getNbrTRNA(inputfile):
    """Gets the number of all tRNA present in the input file.

    :param inputfile: name of the file containing the list of informations
        for each tRNA.
    :type inputfile: string

    :returns: listepiRNA, liste of id piRNA.
    :rtype: list
    """
    listeTRNA = []
    inputfile = open(inputfile,"r")
    for line in inputfile:
        if (re.search('^[0-9]', line)):
            words = line.split('\t')
            infoTRNA = words[1]
            gene = infoTRNA.split('(')[0].replace(" ","")
            if gene not in listeTRNA:
                listeTRNA.append(gene)
    inputfile.close()
    return listeTRNA

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'G4AnnotationTRNA')
    parser.add_argument ('-i', '--inputfile', default = \
        '/home/anais/Documents/Data/Human/G4Conserve-master/' + \
        'tRNA/hg19-mature-tRNAs.csv')
    parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
    parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
    parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
    parser.add_argument ('-E', '--EXTENSION', default = 100)
    parser.add_argument ('-W', '--WINDOW', default = 60)
    parser.add_argument ('-S', '--STEP', default = 10)
    return parser

def main(inputfile, dicoParam):
    G4DetectedIntRNA = ReturnG4InGene(inputfile, dicoParam)
    length = getLengthTRNA(inputfile)
    tRNAWithPG4 = getNbrTRNAWithPG4(G4DetectedIntRNA)
    tRNA = getNbrTRNA(inputfile)

    print '\nInfo PG4 detected in tRNA:'
    print '\tName File  :'+inputfile
    if G4DetectedIntRNA:
        for key, values in G4DetectedIntRNA.items():
            print key+'\t'+'\t'.join(values)
    else:
        print '        NO PG4 DETECTED'
    print '\nNumber PG4 detected in tRNA:    ' + str(len(G4DetectedIntRNA))
    print 'Number tRNA with PG4:    ' + str(len(tRNAWithPG4))
    print 'Number tRNA tested:    ' + str(len(tRNA))
    print 'Length of all tRNA (nt):    ' + str(length)

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    inputfile = arg.inputfile
    dicoParam = rF.createDicoParam(arg)
    main(inputfile, dicoParam)
