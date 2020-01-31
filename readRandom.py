#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""
``init.py`` **module description**:
This module has as input the output of G4RNAScreener and a dictionnary
with all the parameters that were used for G4RNA Screener.
This script filter the output of G4RNA Screener to only keep windows
over thresholds. Then all those windows are merged if they are
overlapping. Overlapping windows are pG4.
.. moduleauthor:: Anaìs Vannutelli, Michel Garant, Sarah Bellamiti, Jean-Pierre Perreault and Aida Ouangraoua
March 2019
Université de Sherbrooke Canada
Laboratoty CoBiUS and Jean-Pierre Perreault
"""

import os
import re
import argparse
import pandas as pd
from pprint import pprint
import recurentFunction as rF

def writepG4Shuffle(dico, path, chr):
	output = open(path+'/pG4r_shuffle_'+chr+'.csv','w')
	output.write('pG4rID\tcGcC\tG4Hunter\tsequenceG4\tG4NN\n')
	for pG4r in dico :
		output.write(pG4r+'\t'+'\t'.join(dico[pG4r])+'\n')

def readLineG4ScreenerLoc(line):
	""" Transcforms a line from the output of G4RNA screener into a dictionary.

	There are many informations in the output of G4RNA screener so
	they will be contains in a dictionary for more readability.

	:paran line: one line from the output G4RNA screener.
	:type line: string

	:returns: dicoLine, line of G4RNA Screener parsed.
	:rtype: dictionary
	"""
	words = line.split('\t') # get all element in a list
	if words[7] == '\n':
		words[7] = '0.0'
	dicoLine = {'Description' : words[1],
				'GeneID' : words[1].split(':')[0],
				'Strand' : words[1].split(':')[4],
				'cGcC' : float(words[2]),
				'g4H' : float(words[3]),
				'WindowSeq' : words[4],
				'WindowStart' : int(words[5]),
				'WindowEnd' : int(words[6]),
				'g4NN' : float(words[7]),
				'locationStart' : int(words[1].split(':')[3].split('-')[0]),
				'locationEnd' : int(words[1].split(':')[3].split('-')[1])}
	return dicoLine

def readLineG4ScreenerJun(line, feature):
	""" Transcforms a line from the output of G4RNA screener into a dictionary.

	There are many informations in the output of G4RNA screener so
	they will be contains in a dictionary for more readability.

	:paran line: one line from the output G4RNA screener.
	:type line: string
	:param StrandByGene: contain al strand of all genes.
	:type StrandByGene: dictionary
	:param feature: name of the feature (gene or junction).
	:type feature: string

	:returns: dicoLine, line of G4RNA Screener parsed.
	:rtype: dictionary
	"""
	words = line.split('\t') # get all element in a list
	dicoLine = {'Description' : words[1],
				'GeneID' : words[1].split(':')[0],
				'Strand' : words[1].split(':')[4],
				'cGcC' : float(words[2]),
				'g4H' : float(words[3]),
				'WindowSeq' : words[4],
				'WindowStart' : int(words[5]),
				'WindowEnd' : int(words[6]),
				'g4NN' : float(words[7]),
				'locationStart' : int(words[1].split(':')[3].split('-')[0]),
				'locationEnd' : int(words[1].split(':')[3].split('-')[1])}
	return dicoLine

def getG4(G4DetectedInGene, inputfile, dicoParam):
	"""Merges all windows from a gene upper thresholds, it will be predicted G4.

	This function browse all windows returned by G4RNA Screener and
	keep only those over the thresholds. It also merge overlapping
	windows. If two windows upper thresholds are separated with one window
	which is under thresholds, the two windows will not be merge and will be
	concidered as 2 pG4.

	:param G4DetectedInGene: pG4 with its scores and sequence, predicted in
		genes.
	:type G4DetectedInGene: dictionary
	:param inputfile: name of the outputfile of G4RNA screener.
	:type inputfile: string
	:param dicoParam: all parameters given to G4RNA screener.
	:type dicoParam: dictionary

	:returns: G4DetectedInGene, contains all predicted G4 in genes,
		with there scores.
	:rtype: dictionary
	"""
	oldPassed = False # boolean about the previous windows
	passed = False # true if over threshold, elsewise false
	inputfile = open(inputfile,'r')
	for line in inputfile:
		if (re.search('^[0-9]', line)):
			# if the line is not the header of the file
			dicoLine = readLineG4ScreenerLoc(line)
			if (dicoLine['cGcC'] >= dicoParam['cGcC'] and
				dicoLine['g4H'] >= dicoParam['g4H'] and
				dicoLine['g4NN'] >= dicoParam['g4NN']):
				# if window over thresholds
				passed = True
				if (oldPassed != passed):
					# first windows over the threshold
					descriptionOverThreshold = dicoLine['Description']
					# assignation of description for this window
					sequenceG4 = dicoLine['WindowSeq']
					oldPassed = passed
					if dicoLine['Strand'] == '1':
						startG4 = dicoLine['WindowStart']
						endG4 = dicoLine['WindowEnd']
					else:
						startG4 = dicoLine['locationEnd'] - \
								(dicoLine['WindowStart'] -\
								dicoLine['locationStart'])
						endG4 = startG4 - (len(dicoLine['WindowSeq'])) + 1
					listeCGcC = [ dicoLine['cGcC'] ]
					listeG4Hunter = [ dicoLine['g4H'] ]
					listeG4NN = [ dicoLine['g4NN'] ]
				else:
					# not the first windows above the thresholds
					sequenceG4 = rF.addWindowToG4Seq(sequenceG4,
								dicoLine['WindowSeq'],
								dicoParam['Step'], dicoParam['Window'])
					if dicoLine['Strand'] == '1' :
						endG4 = dicoLine['WindowEnd']
					else :
						endG4 = startG4 - (len(sequenceG4)) + 1
					listeCGcC.append(dicoLine['cGcC'])
					listeG4Hunter.append(dicoLine['g4H'])
					listeG4NN.append(dicoLine['g4NN'])
			if (dicoLine['cGcC'] < dicoParam['cGcC'] or
				dicoLine['g4H'] < dicoParam['g4H'] or
				dicoLine['g4NN'] < dicoParam['g4NN'] or
				descriptionOverThreshold != dicoLine['Description']):
				# one of the score is under his threshold
				# or this windows is from another gene
				passed = False # update
				if (oldPassed != passed ):
					# last windows before under the thresolds
					meanCGcC = rF.mean(listeCGcC)
					meanG4Hunter = rF.mean(listeG4Hunter)
					meanG4NN = rF.mean(listeG4NN)
					oldPassed = passed # update
					if descriptionOverThreshold != dicoLine['Description']:
						headerG4 = rF.createIdpG4rShuffle(descriptionOverThreshold,
									startG4, endG4, dicoLine['Strand'])
					else:
						headerG4 = rF.createIdpG4rShuffle(dicoLine['Description'],
									startG4, endG4, dicoLine['Strand'])
					if headerG4 not in G4DetectedInGene and dicoLine['Strand']:
						G4DetectedInGene[headerG4] = [str(meanCGcC),
													str(meanG4Hunter),
													sequenceG4,
													str(meanG4NN)]
	inputfile.close()
	return G4DetectedInGene

def getChromosomalPositionForJunction(coord, strand, junLength,
									startIntron, endIntron):
	"""Computes chromosomal coordinates for a.

	This function will leads to a function depending on the strand of the
	junction.

	:param coord: start of end of a junction.
	:type coord: integer
	:param strand: strand of the junctino (1 or -1).
	:type strand: string
	:param junLength: length of the junction from each side (by default 100).
	:type junLength: integer
	:param startIntron: start of the intron from the junction.
	:type startIntron: integer
	:param endIntron: end of the intron from the junction.
	:type endIntron: integer

	:returns: coord, the new coordinates.
	:rtype: integer
	"""
	if strand == '1':
		coord = getChromosomalPositionForwardStrand(coord,
				junLength, startIntron, endIntron)
	elif strand == '-1':
		coord = getChromosomalPositionReverseStrand(coord,
				junLength, startIntron, endIntron)
	else :
		coord = None
	return coord

def isG4OnJunction(startG4, endG4, startBorder, endBorder):
	"""Determines if the G4 is on a junction.

	This function search if there is an overlap between the pG4 and sequences
	upstream and downstream of an intron. If the pG4 overlaps both, then it is
	on a junction.

	:param startG4: start of a pG4.
	:type StartG4: integer
	:param endG4: end of a pG4.
	:type endG4: integer
	:startBorder: start of an intron
	:type startBorder: integer
	:param endBorder: end of an intron
	:type endBorder: integer

	:returns: onJunction, true if the pG4 is on a junction, else false.
	:rtype: boolean
	"""
	if ((startG4 < startBorder and endG4 > startBorder) or
		(startG4 > startBorder and endG4 < startBorder)):
		onJunction = True
	else:
		onJunction = False
	return onJunction

def getChromosomalPositionForwardStrand(position, EXTENSION,
									startIntron, endIntron):
	"""Computes the chromosomal coordinatea from a forward junction.

	Gives real position (chromosomal) of start or end from a G4
	around a junction.

	:param position: start of end of a junction.
	:type position: integer
	:param EXTENSION: length of the junction from each side (by default 100).
	:type EXTENSION: integer
	:param startIntron: start of the intron from the junction.
	:type startIntron: integer
	:param endIntron: end of the intron from the junction.
	:type endIntron: integer

	:returns: position, the new coordinates.
	:rtype: integer
	"""
	if position <= EXTENSION:
		# upstream squence (before the junction)
		position = startIntron - 1 - EXTENSION + position
	else:
		# downstream sequence (after the junction)
		position = endIntron + 1 - EXTENSION + position - 1
	return position

def getChromosomalPositionReverseStrand(position, EXTENSION,
									startIntron, endIntron):
	"""Computes the chromosomal coordinate from a reverse junction.

	Give real position (chromosomal) of start or end from a G4
	around a junction.

	:param position: start of end of a junction.
	:type position: integer
	:param EXTENSION: length of the junction from each side (by default 100).
	:type EXTENSION: integer
	:param startIntron: start of the intron from the junction.
	:type startIntron: integer
	:param endIntron: end of the intron from the junction.
	:type endIntron: integer

	:returns: position, the new coordinate.
	:rtype: integer
	"""
	if position <= EXTENSION:
		# upstream squence (before the junction)
		position = startIntron + 1 + EXTENSION + 1 - position
	else:
		# downstream sequence (after the junction)
		position = endIntron + EXTENSION - position
	return position

def returnG4InJunction(G4DetectedInJunction, inputfile, dicoParam):
	"""Merges all windows from junctions upper thresholds, it will be pG4.

	This function browses all windows returned by G4RNA Screener and
	keep only those over the thresholds. It also merge overlapping
	windows. If two windows upper thresholds are separated with one window
	which is under thresholds, the two windows will not be merge and will be
	concidered as 2 pG4.

	:param G4DetectedInJunction: pG4 with its scores and sequence, those pG4
		are predicted in junctions.
	:type G4DetectedInJunction: dictionary
	:param inputfile: name of the outputfile of G4RNA screener.
	:type inputfile: string
	:param dicoParam: all parameters given to G4RNA screener.
	:type dicoParam: dictionary

	:returns: G4DetectedInJunction, contains all pG4 in juctions,
		with there scores.
	:rtype: dictionary
	"""
	oldPassed = False
	passed = False
	descriptionOverThreshold = ''
	inputfile = open(inputfile,'r')
	for line in inputfile:
		if (re.search('^[0-9]', line)):
			# if the line is not the header of the file
			dicoLine = readLineG4ScreenerJun(line, 'Junction')
			if (dicoLine['cGcC'] >= dicoParam['cGcC'] and
				dicoLine['g4H'] >= dicoParam['g4H'] and
				dicoLine['g4NN'] >= dicoParam['g4NN'] and
				dicoLine['WindowStart']):
				# window over thresholds
				onJunction = False
				passed = True
				if (oldPassed != passed):
					# first windows over thresholds
					descriptionOverThreshold = dicoLine['Description']
					sequenceG4 = dicoLine['WindowSeq']
					oldPassed = passed # update
					startFirstWindow = dicoLine['locationStart']
					endFirstWindow = dicoLine['locationEnd']
					startG4 = dicoLine['WindowStart']
					endG4 = dicoLine['WindowEnd']
					listeCGcC = [ dicoLine['cGcC'] ]
					listeG4Hunter = [ dicoLine['g4H'] ]
					listeG4NN = [ dicoLine['g4NN'] ]
				else:
					# not the first windows above the thresholds
					sequenceG4 = rF.addWindowToG4Seq(sequenceG4,
								dicoLine['WindowSeq'],
								dicoParam['Step'],
								dicoParam['Window'])
					endG4 = dicoLine['WindowEnd']
					listeCGcC.append(dicoLine['cGcC'])
					listeG4Hunter.append(dicoLine['g4H'])
					listeG4NN.append(dicoLine['g4NN'])
			if (dicoLine['cGcC'] < dicoParam['cGcC'] or
				dicoLine['g4H'] < dicoParam['g4H'] or
				dicoLine['g4NN'] < dicoParam['g4NN'] or
				descriptionOverThreshold != dicoLine['Description']):
				# one of the score is under his threshold
				# or this windows is from another gene
				passed = False
				if (oldPassed != passed ):
					# last windows before under the thresolds
					meanCGcC = rF.mean(listeCGcC)
					meanG4Hunter = rF.mean(listeG4Hunter)
					meanG4NN = rF.mean(listeG4NN)
					oldPassed = passed
					if (startG4 > 40 or endG4 < 160):
						if descriptionOverThreshold != dicoLine['Description']:
							description = descriptionOverThreshold.split(':')[0]+':junction:'+':'.join(descriptionOverThreshold.split(':')[2:])
						else:
							description = dicoLine['Description'].split(':')[0]+':junction:'+':'.join(dicoLine['Description'].split(':')[2:])
						if description not in G4DetectedInJunction:
							G4DetectedInJunction[description] = [str(meanCGcC),
															str(meanG4Hunter),
															sequenceG4,
															str(meanG4NN)]
	inputfile.close()
	return G4DetectedInJunction

def main(dicoParam, path, chr):
	G4DetectedInGene = {}
	G4DetectedInJunction = {}
	directoryCSV = path + '/Data/'+chr+'/CSVFile'
	# get g4 from the ouput of G4RNA Screener
	for path2, dirs, files in os.walk(directoryCSV):
		# for each element of the directory to passed
		for filename in files: # for each files (all csv)
			inputfile = directoryCSV + '/' + filename
			if ('location' in filename):
				G4DetectedInGene = getG4(G4DetectedInGene,
									inputfile, dicoParam)
			elif ('junction' in filename):
				G4DetectedInJunction = returnG4InJunction(G4DetectedInJunction,
										inputfile, dicoParam)
	pG4inRegion = G4DetectedInGene
	pG4inRegion.update(G4DetectedInJunction)
	writepG4Shuffle(pG4inRegion, path+'/Results/perChromosome', chr)

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'readRandom')
	parser.add_argument ('-p', '--path', default = '/home/vana2406/scratch/TestRepro')
	parser.add_argument ('-chr', '--Chromosome', default = 'Y')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	chr = arg.Chromosome
	path = arg.path
	path = path
	dicoParam = rF.createDicoParam(arg)
	main(dicoParam, path, chr)
	print("\tDone")
