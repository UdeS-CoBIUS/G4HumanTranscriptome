#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	August 2018

Description:
	This script aims to predict pG4r in circular RNA and to computes
		statistiques about it. Firstly, windows from G4RNA screener are merge
		into pG4r. Then density is computed. The output is printed on the
		screen.

Data availability:
	* fasta sequences used with G4RNA screener were downloaded from circBase in August 2018.
	* G4RNA screener was launch in August 2018. Results may change if G4RNA is update later.

"""

import sys
import os
import re
import argparse
import recurentFunction as rF

def readLineCircu(words):
	"""Reads a line from circBase file.

	:param words: contains all informations from a line of circBase.
	:type words: list

	:returns: dicoL, contains all informations of circBase but parsed.
	:rtype: dictionary
	"""
	dicoL = {'circuInfo' : words[1],
			'cGcC' : float(words[2]),
			'g4H' : float(words[3]),
			'Sequence' : words[4],
			'wStart' : int(words[5]),
			'wEnd' : int(words[6]),
			'g4NN' : float(words[7].rstrip())}
	infoSplit = dicoL['circuInfo'].split('|')
	positions = infoSplit[1].split(':')
	dicoL.update({'circuId' : infoSplit[0],
				'Chromosome' : positions[0],
				'circuStart' : int(positions[1].split('-')[0]),
				'circuEnd' : int(positions[1].split('-')[1].replace("+","")),
				'Strand' : positions[1][-1]})
	dicoL['Description'] = dicoL['Chromosome'] + '|' + str(dicoL['circuStart'])\
		+ '-' + str(dicoL['circuEnd'])
	return dicoL

def ReturnG4InCircu(G4Detected, inputfile, dicoParam):
	"""Merges all windows from a gene upper thresholds, it will be predicted G4.

	This function browses all windows returned by G4RNA Screener and
	keep only those over the thresholds. It also merges overlapping
	windows. If two windows upper thresholds are separated with one window
	which is under thresholds, the two windows will not be merge and will be
	concidered as 2 pG4.

	:param G4Detected: pG4r with its scores and sequence, predicted in
		genes.
	:type G4Detected: dictionary
	:param inputfile: name of the outputfile containing the liste of
		informations for each circular.
	:type inputfile: string
	:param dicoParam: all parameters given to G4RNA screener.
	:type dicoParam: dictionary

	:returns: G4Detected, contains all predicted G4 in genes, with there scores.
	:rtype: dictionary
	"""
	rnaType = 'circu'
	location = 'ExonNC'
	oldPassed = False
	passed = False
	overtakingSequence = False
	first = False
	tmp = []
	upstream = ''
	downstream = ''
	condition = False
	descriptionOverThreshold = ''
	inputfile = open(inputfile,"r")
	for line in inputfile:
		if re.search('^[0-9]', line):
			words = line.split('\t')
			if len(words) == 8:
				dL = readLineCircu(words)
				if (dL['cGcC'] >= dicoParam['cGcC'] and
					dL['g4H'] >= dicoParam['g4H'] and
					dL['g4NN'] >= dicoParam['g4NN']):
					passed = True
					if (oldPassed != passed):
						# if it's the first windows, beginning if the passage
						descriptionOverThreshold = dL['Description']
						sequenceG4 = dL['Sequence']
						oldPassed = passed
						if dL['Strand'] == '+':
							startG4 = dL['wStart']
							# same as startWindow of G4 screener
							endG4 = dL['wEnd']
							# same as endWindow of G4 screener
						else:
							startG4 = dL['circuEnd'] - dL['wStart'] + 1
							endG4 = startG4 - len(dL['Sequence']) + 1
						listeCGcC = [ dL['cGcC'] ]
						listeG4Hunter = [ dL['g4H'] ]
						listeG4NN = [ dL['g4NN'] ]
					else:
						# if it's not the first windows above the thresholds
						sequenceG4 = rF.addWindowToG4Seq(sequenceG4,
								dL['Sequence'], dicoParam['Step'],
								dicoParam['Window'])
						if dL['Strand'] == '+':
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
					if (oldPassed != passed ):
						# last windows before under the thresolds
						meanCGcC = rF.mean(listeCGcC)
						meanG4Hunter = rF.mean(listeG4Hunter)
						meanG4NN = rF.mean(listeG4NN)
						oldPassed = passed # update
						headerG4 = rF.createIdG4(dL['circuId'],
								startG4, endG4, dL['Strand'])
						if dL['Strand'] == '+':
							if (startG4 > dL['circuStart'] and
								startG4 < dL['circuEnd'] and
								endG4 > dL['circuEnd']):
								# overlap
								condition = True
						else:
							if (startG4 > dL['circuStart'] and
								startG4 < dL['circuEnd'] and
								endG4 < dL['circuStart']): # overlap
								condition = True
						if condition == True: # if overlap
							location = 'junction'
							if dL['Strand'] == '+':
								endG4 = dL['circuStart'] + \
									(endG4 - dL['circuEnd'])
							else:
								endG4 = dL['circuEnd'] - \
									(dL['circuStart'] - endG4)
							headerG4 = rF.createIdG4(dL['circuId'],
									startG4, endG4, dL['Strand'])
							G4Detected[headerG4] = [str(meanCGcC),
													str(meanG4Hunter),
													sequenceG4, str(meanG4NN),
													location, rnaType]
						else:
							if (G4Detected.has_key(headerG4) == False and
								dL['Strand'] != None):
								G4Detected[headerG4] = [str(meanCGcC),
														str(meanG4Hunter),
														sequenceG4,
														str(meanG4NN),
														location, rnaType]
						condition = False
						location = 'ExonNC'
	inputfile.close()
	return G4Detected


def ReturnLengthInCircu(Length, inputfile):
	"""Computes the length of all circRNA.

	:param Length: length of each circular regions.
	:type Length: dictionary
	:param inputfile: name of files containing the list of informations
	 	for each circular RNA.
	:type inputfile: string

	:returns: Length, updated.
	:rtype: dictionary
	"""
	inputfile = open(inputfile,"r")
	for line in inputfile:
		if (re.search('^[0-9]', line)):
			words = line.split('\t')
			if (len(words) == 8):
				numRow = words[0]
				infoCircu = words[1]
				idCircu = infoCircu.split('|')[0]
				startBorder = \
					infoCircu.split('|')[1].split(':')[1].split('-')[0]
				endBorder = \
					infoCircu.split('|')[1].split(':')[1].split('-')[1]\
					.replace("+","")
				if idCircu not in Length:
					Length[idCircu] = \
						rF.GetLengthFraction(startBorder, endBorder)
	inputfile.close()
	return Length

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'G4AnnotationCirculaire')
	parser.add_argument ('-specie', '--specie', default = 'HS')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

def main(specie, dicoParam):
	directory = '/home/anais/Documents/Data/Human/G4Conserve-master/' + \
		'circulaire/Circu/'
	G4DetectedInCircu = {}
	LengthInCircu = {}
	length = 0
	for path, dirs, files in os.walk(directory):
		for filename in files:
			inputfile = directory + filename
			G4DetectedInCircu = ReturnG4InCircu(G4DetectedInCircu, inputfile,
								dicoParam)
			LengthInCircu = ReturnLengthInCircu(LengthInCircu, inputfile)
	print '----------------------------------------'
	sommeLength = 0
	for key, values in LengthInCircu.items():
		sommeLength = sommeLength + values
	print sommeLength , len(G4DetectedInCircu)
	print len(G4DetectedInCircu) / float(sommeLength) * 1000
	print '---------------------------------------------'
	print '\npG4r number detected in circRNA:	' + str(len(G4DetectedInCircu))
	print 'Total number of circRNA tested:	' + str(len(LengthInCircu))
	print 'Length of all circRNA (nt):	' + str(sommeLength)
	print 'Density in Kb :	' + str(round(len(G4DetectedInCircu)/ \
                                float(sommeLength) * 1000, 2))
	print 'Percent circRNA with Kb:	' + str(round(len(G4DetectedInCircu)/ \
                                        float(len(LengthInCircu)) * 100,2))

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	specie = arg.specie
	dicoParam = rF.createDicoParam(arg)
	main(specie, dicoParam)
