#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	December 2017

Description:
	This script aims to annote the pG4 region of one chromosome of one specie.
	From csv files which contain windows from G4RNA screener, this module
	annote the G4 region for each transcript.

Data availability:
	* fasta sequences used with G4RNA screener were download from Ensembl. The human assembly is GRCh38p12.
	* G4RNA screener was launch in August 2018. Results may change if G4RNA is update later.
	* biotypes (classes and subclasses) are downloaded from Ensembl via biomart in August 2018.
	* transcripts informations are also downloaded from Ensembl in August 2018	and parsed by ReplaceInformation.py.

"""

import sys
import os
import re
import recurentFunction as rF
from pprint import pprint
import argparse

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

def readLineG4Screener(line, StrandByGene, feature):
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
				'GeneID' : words[1].split('|')[0],
				'Strand' : StrandByGene.get(words[1].split('|')[0]),
				'cGcC' : float(words[2]),
				'g4H' : float(words[3]),
				'WindowSeq' : words[4],
				'WindowStart' : int(words[5]),
				'WindowEnd' : int(words[6]),
				'g4NN' : float(words[7])}
	if feature == 'Gene':
		dicoLine.update({feature+'Start' : int(words[1].split('|')[1]),
						feature+'End' : int(words[1].split('|')[2])})
	else:
		# for a junction, it is the coord of the intron
		dicoLine.update({feature+'Start' : int(words[1].split('|')[1]\
						.split('-')[0]),
						feature+'End' : int(words[1].split('|')[1]\
						.split('-')[1])})
	return dicoLine

def addCoordInlist(coordFeatureList, storageList):
	"""Stores coordinate from a list of feature into an other list.

	Browses the initial list of coords of a feature, parse them and
	then add them in a list. This list stores all coords of a transcript (exon,
	intron and UTR).

	:param coordFeatureList: contains coordinates for a feature (intron or exon).
	:type coordFeatureList: list
	:param storageList: where the parsed coord are stored.
	:type storageList: list

	:returns: storageList, list of coordinates updated.
	:rtype: list
	"""
	for coord_feature in coordFeatureList:
		if coord_feature:
			startFeature = coord_feature.split('-')[0]
			endFeature = coord_feature.split('-')[1]
			storageList.append(int(startFeature))
			storageList.append(int(endFeature))
	return storageList

def getCoordOfTranscript(start5 , end5 , start3 , end3,
							exonList, intronList, strand):
	"""Retrieves start and end of a transcript.

	:param start5: start of the 5'UTR.
	:type start5: integer
	:param end5: end of the 5'UTR.
	:type end5: integer
	:param start3: start of the 3'UTR.
	:type start3: integer
	:param end3: end of the 3'UTR.
	:type end3: integer
	:param exonList: contain coordinates of all exons from a transcript.
	:type exonList: list
	:param intronList: contain coordinates of all introns from a transcript.
	:type intronList: list
	:param strand: strand of a transcript.
	:type strand: string

	:returns: [inferiorCoord, superiorCoord], start and end of a transcript.
	:rtype: list
	"""
	listCoords = [] # store coords the transcript's features
	if not start5:
		# if the transcript doesn't have a 5' UTR
		if not end3:
			# if the transcript doesn't have a 3' UTR
			# coords are given by exons
			listCoords = addCoordInlist(exonList, listCoords)
			listCoords = addCoordInlist(intronList, listCoords)
			if strand == '1':
				# if forward gene then 5' < 3'
				inferiorCoord = min(listCoords)
				superiorCoord = max(listCoords)
			else:
				# if reverse gene then 5' > 3'
				inferiorCoord = max(listCoords)
				superiorCoord = min(listCoords)
		else:
			# if the transcript doesn't have 5' UTR
			# coords are given by exons and 3'
			listCoords = addCoordInlist(exonList, listCoords)
			listCoords = addCoordInlist(intronList, listCoords)
			if strand == '1':
				inferiorCoord = min(listCoords)
			else:
				inferiorCoord = max(listCoords)
			superiorCoord = int(end3)
	elif not start3: # if the transcript doesn't have a 3'UTR
		# border will be given by exons and 5'
		listCoords.append(start5)
		listCoords.append(end5)
		listCoords = addCoordInlist(exonList, listCoords)
		listCoords = addCoordInlist(intronList, listCoords)
		if strand == '1':
			superiorCoord = max(listCoords)
		else:
			superiorCoord = min(listCoords)
		inferiorCoord=int(start5)
	else: # if the transcript has a 5'UTR, a 3'UTR, exons and introns
		inferiorCoord = int(start5)
		superiorCoord = int(end3)
	return [inferiorCoord, superiorCoord]

def isG4InTranscript(strand, coordG4, coordsTranscript):
	"""Searches if the pG4 is on the transcript.

	:param strand: strand of the transcript and of the pG4.
	:type strand: string
	:param coordG4: start and end of a pG4.
	:type coord: list
	:param coordsTranscript: start and end of a transcript.
	:type coordsTranscript: list

	:returns: inTranscript, true if the pG4 is on the transcript else false.
	:rtype: boolean
	"""
	startG4 = int(coordG4[0])
	endG4 = int(coordG4[1])
	inferiorCoord = int(coordsTranscript[0])
	superiorCoord = int(coordsTranscript[1])
	inTranscript = False
	if strand == '1': # forward strand
		if (startG4 >= inferiorCoord and
			startG4 <= superiorCoord and
			endG4 >= inferiorCoord and
			endG4 <= superiorCoord):
			inTranscript = True
	else: # reverse strand
		if (startG4 <= inferiorCoord and
			startG4 >= superiorCoord and
			endG4 <= inferiorCoord and
			endG4 >= superiorCoord):
			inTranscript = True
	return inTranscript

def getLocationFeatureReverse(feature, featureList, startG4, endG4, biotype):
	"""Retrieves the location of the feature on the reverse strand.

	:param feature: name of the feature (gene or junction).
	:type feature: string
	:param featureList: all coords of featured elements.
	:type featureList: list
	:param startG4: start of the pG4.
	:type startG4: integer
	:param endG4: end of the pG4.
	:type endG4: integer
	:param biotype: biotype of the transcript.
	:type biotype: string

	:returns: location, the location of the pG4.
	:rtype: string
	"""
	for coord_feature in featureList:
			if coord_feature:
				# if transcript contain exon
				featureStart = int(coord_feature.split('-')[0])
				featureEnd = int(coord_feature.split('-')[1])
				if endG4 <= featureStart and endG4 >= featureEnd:
					# pG4 on an intron
					if startG4 >= featureEnd:
						# pG4 entierly in an intron
						location = biotype
						if feature == 'Intron' and biotype == 'CDS':
							location = 'Intron'
						elif feature == 'Intron and ' and biotype != 'CDS':
							location = 'IntronNC'
					else:
						# pG4 overlapping an exon and an intron
						if biotype == 'CDS':
							if feature == 'Exon':
								location = 'junction_'+biotype+'_Intron'
							else:
								location = 'junction_Intron_'+biotype
						else:
							if feature == 'Exon':
								location = 'junction_'+biotype+'_IntronNC'
							else:
								location = 'junction_IntronNC_'+biotype
	return location

def getLocationFeatureForward(feature, featureList, startG4, endG4, biotype):
	"""Retrieves the location of a feature on the forward strand.

	:param feature: name of the feature (gene or junction).
	:type feature: string
	:param featureList: all coords of featured elements for this transcript.
	:type featureList: list
	:param startG4: start of the pG4.
	:type startG4: integer
	:param endG4: end of the pG4.
	:type endG4: integer
	:param biotype: biotype of the transcript.
	:type biotype: string

	:returns: location, the location of a pG4.
	:rtype: string
	"""
	for coord_feature in featureList:
			if coord_feature:
				# if transcript contain exon
				featureStart = int(coord_feature.split('-')[0])
				featureEnd = int(coord_feature.split('-')[1])
				if startG4 >= featureStart and startG4 <= featureEnd:
					# pG4 on an intron
					if endG4 <= featureEnd:
						# pG4 entierly in an intron
						location = biotype
						if feature == 'Intron' and biotype == 'CDS':
							location = 'Intron'
						elif feature == 'Intron and ' and biotype != 'CDS':
							location = 'IntronNC'
					else:
						# pG4 overlapping an exon and an intron
						if biotype == 'CDS':
							if feature == 'Exon':
								location = 'junction_'+biotype+'_Intron'
							else:
								location = 'junction_Intron_'+biotype
						else:
							if feature == 'Exon':
								location = 'junction_'+biotype+'_IntronNC'
							else:
								location = 'junction_IntronNC_'+biotype
	return location

def getlocationForward(BiotypeByTranscript, ProteinCoding, transcriptId,
						coordG4, end5, start3, exonList, intronList):
	"""Gets the locaion of pG4 region in a transcript on forward strand.

	:param BiotypeByTranscript: contain all bitoype of all transcripts.
	:type BiotypeByTranscript: dictionary
	:param ProteinCoding: contain all biotype names of protein coding.
	:type ProteinCoding: list
	:param transcriptId: ensembl id of the transcript.
	:type transcriptId: string
	:param coordG4: coords of the pG4
	:type coordG4: list
	:param end5: end of the 3'UTR.
	:type end5: integer
	:param start3: start of the 5'UTR.
	:type start3: integer
	:param exonList: all coords of exons for this transcript.
	:type exonList: list
	:param intronList: all coords of introns for this transcript.
	:type intronList: list

	:returns: location, the location of the pG4.
	:rtype: string
	"""
	startG4 = int(coordG4[0])
	endG4 = int(coordG4[1])
	location = ''
	biotypeTranscript = BiotypeByTranscript.get(transcriptId)
	if biotypeTranscript in ProteinCoding:
		codant = 'CDS'
	else:
		codant = 'ExonNC'
	if end5 and startG4 <= int(end5): # if pG4 is positioned on the 5'UTR
		if endG4 <= int(end5): # if pG4 is enterly in the 5'UTR
			location='5'
		else: # if G4 is overlapping the 5'UTR
			location = 'junction_5_'+codant
	elif start3 and endG4 >= int(start3): # if pG4 is positioned on the 3'UTR
		if startG4 >= int(start3): # if pG4 is enterly in the 3'UTR
			location = '3'
		else: # if pG4 is overlapping the 3'UTR
			location = 'junction_'+codant+'_3'
	else: # if pG4 is positioned in an intron
		location = getLocationFeatureForward('Exon', exonList,
						startG4, endG4, codant)
		location = getLocationFeatureForward('Intron', intronList,
						startG4, endG4, codant)
	return location

def getlocationReverse(BiotypeByTranscript,
						ProteinCoding,
						transcriptId,
						coordG4, end5, start3,
						exonList, intronList):
	"""Returns the location of pG4 region in a transcript on reverse strand.

	:param BiotypeByTranscript: contain all bitoype of all transcripts.
	:type BiotypeByTranscript: dictionary
	:param ProteinCoding: contain all biotype names of protein coding.
	:type ProteinCoding: list
	:param transcriptId: ensembl id of the transcript.
	:type transcriptId: string
	:param coordG4: coords of the pG4
	:type coordG4: list
	:param end5: end of the 3'UTR.
	:type end5: integer
	:param start3: start of the 5'UTR.
	:type start3: integer
	:param exonList: all coords of exons for this transcript.
	:type exonList: list
	:param intronList: all coords of introns for this transcript.
	:type intronList: list

	:returns: location, the location of the pG4.
	:rtype: string
	"""
	startG4 = int(coordG4[0])
	endG4 = int(coordG4[1])
	location = 'NAN'
	biotypeTranscript = BiotypeByTranscript.get(transcriptId)
	if biotypeTranscript in ProteinCoding:
		codant = 'CDS'
	else:
		codant = 'ExonNC'
	if end5 and endG4 >= int(end5): # if pG4 is positioned on the 5'UTR
		if startG4 >= int(end5): # if pG4 is enterly in the 5'UTR
			location = '5'
		else:	# if G4 is positioned in the 5' region with overlapping
			location = 'junction_5_'+codant
	elif start3 and startG4 <= int(start3): # if pG4 is positioned on the 3'UTR
		if endG4 <= int(start3): # if pG4 is enterly in the 3'UTR
			location = '3'
		else: # if pG4 is overlapping the 3'UTR
			location = 'junction_'+codant+'_3'
	else: # if pG4 is positioned in an intron
		location = getLocationFeatureReverse('Exon', exonList,
						startG4, endG4, codant)
		location = getLocationFeatureReverse('Intron', intronList,
						startG4, endG4, codant)
	return location

def getlocationG4InJunction(BiotypeByTranscript,
								ProteinCoding,
								transcriptId):
	"""Sets the location name depending on the biotype.

	If the biotype is one of the codent it will be a CDS, else wise it will be
	ExonNC. This done for a junction.

	:param BiotypeByTranscript: contain biotype of all transcripts.
	:type BiotypeByTranscript: dictionary
	:param ProteinCoding: contain all type of biotype that correspond to a
		protein coding.
	:type ProteinCoding: list
	:param transcriptId: id of the transcript.
	:type transcriptId: string

	:returns: the name of the junction, depending on the transcript biotype.
	:rtype: string
	"""
	location = 'NAJ'
	biotypeTranscript = BiotypeByTranscript.get(transcriptId)
	if biotypeTranscript in ProteinCoding:
		codant = 'CDS'
	else:
		codant = 'ExonNC'
	return 'junction_'+codant+'_'+codant

def getlocationG4InTranscript(BiotypeByTranscript,
								ProteinCoding,
								transcriptId,
								coordG4, end5, start3,
								exonList, intronList, strand):
	""" Returns the location of pG4 region in a transcript.

	Depending on the strand, the location of the pG4 is found.

	:param BiotypeByTranscript: contain all bitoype of all transcripts.
	:type BiotypeByTranscript: dictionary
	:param ProteinCoding: contain all biotype names of protein coding.
	:type ProteinCoding: list
	:param transcriptId: ensembl id of the transcript.
	:type transcriptId: string
	:param coordG4: coords of the pG4
	:type coordG4: list
	:param end5: end of the 3'UTR.
	:type end5: integer
	:param start3: start of the 5'UTR.
	:type start3: integer
	:param exonList: all coords of exons for this transcript.
	:type exonList: list
	:param intronList: all coords of introns for this transcript.
	:type intronList: list
	:param strand: strand of the transcipt.
	:type strand: string

	:returns: location, the location of the pG4.
	:rtype: string
	"""
	if strand == '1': # if gene on a forward strand
		location = getlocationForward(BiotypeByTranscript,
											ProteinCoding,
											transcriptId,
											coordG4, end5, start3,
											exonList, intronList)
	else: # if gene on a reverse strand
		location = getlocationReverse(BiotypeByTranscript,
											ProteinCoding,
											transcriptId,
											coordG4, end5, start3,
											exonList, intronList)
	return location

def addG4InTranscriptome(G4InTranscript, transcriptId,
						descriptionG4, informationsOfG4,
						locationInTranscript, biotypeTranscript):
	"""Stores each pG4 with its information in the transcriptom.

	This fonction add for each pG4 the informations from G4RNA screener
	(cGcC, g4H, sequence and G4NN) to the location of the pG4 and the
	transcript biotype. Those Information are stored in a dictionary.

	:param G4InTranscript: contains all pG4 (key) and all info about it at the
		transcriptom level.
	:type G4InTranscript: dictionary
	:param transcriptId: id of the transcript.
	:type transcriptId: string
	:param descriptionG4: id of the G4 (chromosome, start, end, strand).
	:type descriptionG4: string
	:param informationsOfG4: list of informations obtained by G4-screener
		(cGcC, g4H, sequence and G4NN).
	:type informationsOfG4: list
	:param locationInTranscript: location of G4 on this transcript.
	:type locationInTranscript: string
	:param biotypeTranscript: biotype of this this transcript.
	:type biotypeTranscript: string

	:returns: G4InTranscript, contains id of transcripts with the description of
	 	G4 and their corresponding info of this G4.
	:rtype: dictionary
	"""

	value = []
	for information in informationsOfG4:
		value.append(information)
	value.append(locationInTranscript)
	value.append(biotypeTranscript)
	headerG4 = transcriptId+'|'+descriptionG4 # uniq header pG4
	if headerG4 not in G4InTranscript:
		G4InTranscript[headerG4] = value
	return G4InTranscript

def addG4InGenome(G4InGenome, geneId, descriptionG4, locationInTranscript):
	"""Adds for each gene all pG4 in it and all location possible.

	:param G4InGenome: id of gene with description of G4 and to them correspond
		all location of this G4.
	:type G4InGenome: dictionary
	:param geneId: Ensembl id of the gene.
	:type geneId: string
	:param descriptionG4: id of the G4 (chromosome, start, end, strand).
	:type descriptionG4: string
	:param locationInTranscript: location of G4 on this transcript.
	:type locationInTranscript: string

	:returns: G4InGenome updated.
	:rtype: dictionary
	"""
	headerG4 = geneId+':'+descriptionG4
	if headerG4 not in G4InGenome:
		G4InGenome[headerG4] = locationInTranscript
	else:
		locationInGene=G4InGenome.get(headerG4)
		if locationInGene :
			if locationInTranscript not in locationInGene:
				G4InGenome[headerG4] = locationInGene+';'+locationInTranscript
	return G4InGenome

def addTranscriptPerG4(TranscriptPerG4, descriptionG4, transcriptId):
	"""Adds for a pG4 all transcripts that contain him.

	:param TranscriptPerG4: contains all pG4 and all transcript possessing the
		pG4.
	:type TranscriptPerG4: dictionary
	:param descriptionG4: id of the G4 (chromosome, start, end, strand).
	:type descriptionG4: string
	:param transcriptId: id of the transcript.
	:type transcriptId: string

	:returns: TranscriptPerG4 updated.
	:rtype: dictionary
	"""
	if descriptionG4 not in TranscriptPerG4:
		TranscriptPerG4[descriptionG4] = transcriptId
	else:
		if TranscriptPerG4:
			# dico not empty, he already contains some information
			lastTranscript = TranscriptPerG4.get(descriptionG4)
			TranscriptPerG4[descriptionG4] = lastTranscript +'-'+transcriptId
	return TranscriptPerG4

def isJuncionInTranscript(coordG4,intronList):
	"""Searches if a junction with a pG4 is in a transcript.

	:param coordG4: contains start and end of the pG4.
	:type coordG4: list
	:param intronList: introns contained in this transcript.
	:type intronList: list

	:returns: inTranscript, true if the pG4 is on an intron, false other wise.
	:rtype: boolean
	"""
	startG4 = coordG4[0]
	endG4 = coordG4[1]
	inTranscript = False
	for couple_Intron in intronList:
		if couple_Intron:
			endIntron = int(couple_Intron.split('-')[1])
			if startG4 <= endIntron and endIntron<= endG4:
				inTranscript = True
	return inTranscript

def getUTRCoordByStrand(words, dicoLine):
	"""Retrieves coords of UTR depending on the strand.

	:param words: contain all words of the index file.
	:type words: list
	:param dicoLine: contain all info on a line (except UTR).
	:type dicoLine: dictionary

	:returns: dicoLine, contains all info on a line with UTR.
	:rtype: dictionary
	"""
	if dicoLine['Strand'] == '1':
		dicoLine.update({'start5' : words[7],
						'end5' : words[8],
						'start3' : words[9],
						'end3' : words[10].rstrip()})
	else:
		dicoLine.update({'start5' : words[8],
						'end5' : words[7],
						'start3' : words[10].rstrip(),
						'end3' : words[9]})
	return dicoLine

def readLineIndex(line):
	"""Parses and stores infromations from a line of transcript index.

	The transcript index file contains for each transcript, the gene from where
	it came, its coordinates (start, end, chromosome, strand) and also all
	coordinates of intron/exon/UTR.

	:param line: line from the transcript file.
	:type line: string

	:returns: dicoLine, contains all informations for each transcript.
	:rtype: dictionary
	"""
	words = line.split('|')
	dicoLine = {'idTr' : words[0],
				'idGene' : words[1],
				'Chromosome' : words[2],
				'Strand' : words[3],
				'geneBiotype' : words[4],
				'exonList' : words[5].split(';'),
				'intronList' : words[6].split(';')}
	dicoLine.update(getUTRCoordByStrand(words, dicoLine))
	return dicoLine

def extractionTranscriptPerG4(directory, specie, chromosome, TranscriptPerG4):
	"""Writes in a file all pG4 and all transcript that contains it.

	.. admonition:: Example

		pG4Id1   transcriptId1, transcriptId2, etc

	:param directory: directory for the outpu file.
	:type directory: string
	:param specie: short name of the specie (HS for example).
	:type specie: string
	:param chromosome: chromosome which is analyzed.
	:type chromosome: string
	:param TranscriptPerG4: contains all pG4 and their transcripts
	:type TranscriptPerG4: dictionary
	"""
	output = open(directory + '/' + specie + '_chr' + chromosome + \
				'_TranscriptPerG4.txt', 'w')
	output.write('CoordonnéesG4\tTranscript(s)\n')
	for key, value in TranscriptPerG4.items():
		output.write(key + '\t' + value + '\n')
	output.close()

def extractionG4InGenome(directory, specie, chromosome, G4InGenome):
	"""Writes in a file all pG4 and all transcript that contains it.

	.. admonition:: Example

		geneID|pG4Id   location1, location2, etc

	:param directory: directory for the outpu file.
	:type directory: string
	:param specie: short name of the specie (HS for example).
	:type specie: string
	:param chromosome: chromosome which is analyzed.
	:type chromosome: string
	:param G4InGenome: id of gene with description of G4 and to them correspond
		all location of this G4.
	:type G4InGenome: dictionary
	"""
	output = open(directory+'/'+specie+'_chr'+chromosome+'_G4InGenome.txt','w')
	output.write('InfoG4ByGene\tlocation(s)\n')
	for key,value in G4InGenome.items():
		output.write(key+'\t'+value+'\n')
	output.close()

def extractionG4InTranscript(directory, specie, chromosome, G4InTranscript):
	"""Writes in a file all pG4 and all transcript that contains it.

	.. admonition:: Example

		transcriptId|pG4Id   cGcC G4H sequence G4NN location biotype

	:param directory: directory for the outpu file.
	:type directory: string
	:param specie: short name of the specie (HS for example).
	:type specie: string
	:param chromosome: chromosome which is analyzed.
	:type chromosome: string
	:param G4InTranscript: contains all pG4 (key) and all info about it at the
		transcriptom level.
	:type G4InTranscript: dictionary
	"""
	output = open(directory+'/'+specie+'_chr'+chromosome+'_G4InTranscript.txt','w')
	output.write('InfoG4ByTranscript\tcGcC\tG4Hunter\tsequenceG4\tG4NN\tlocation\ttranscriptBiotype\n')
	for key,value in G4InTranscript.items():
		if None not in value:
			# because some transcriptID from ensembl don't
			# contains info of biotype (as ENST00000604369)
			output.write(key+'\t'+'\t'.join(value)+'\n')
	output.close()

def getInfoAboutpG4(index, BiotypeByTranscript, listeG4InGeneEntire,
					G4InTranscript, G4InGenome, TranscriptPerG4,
					AnnotationTranscript, G4DetectedInGene, ProteinCoding,
					feature):
	"""Searches all locations and biotype of a pG4.

	As a pG4 can be on different transcripts, it can have many location and
	and many biotype. Here, we pass through all pG4 and we find in which
	genes they are, in which transcripts and from what biotype.

	:param index: file name of the index file.
	:type index: string
	:param BiotypeByTranscript: contain all bitoype of all transcripts.
	:type BiotypeByTranscript: dictionary
	:param listeG4InGeneEntire: contain all gene and for each of them a list of
		all the pG4 they have.
	:type listeG4InGeneEntire: dictionary
	:param G4InTranscript: contains all pG4 (key) and all info about it at the
		transcriptom level.
	:type G4InTranscript: dictionary
	:param G4InGenome: id of gene with description of G4 and to them correspond
		all location of this G4.
	:type G4InGenome: dictionary
	:param TranscriptPerG4: contains all pG4 and their transcripts
	:type TranscriptPerG4: dictionary
	:param AnnotationTranscript: information of annotation for all transcripts
		By Ensembl (True if good anotation, else False)
	:type AnnotationTranscript: dictionary
	:param G4DetectedInGene: all pG4 with its score and sequence.
	:type G4DetectedInGene: dictionary
	:param ProteinCoding: contain all biotype names of protein coding.
	:type ProteinCoding: list
	:param feature: name of the feature (junction or gene).
	:type feature: string

	:returns: G4InTranscript, G4InGenome and TranscriptPerG4 are updated.
	:rtype: dictionaries
	"""
	inputfile = open(index,'r') # file opening for reading
	for line in inputfile: # for each transcript
		l = readLineIndex(line)
		biotypeTranscript = BiotypeByTranscript.get(l['idTr'])
		borderTranscript = getCoordOfTranscript(l['start5'], l['end5'],
												l['start3'], l['end3'],
												l['exonList'],
												l['intronList'],
												l['Strand'])
		if l['idGene'] in listeG4InGeneEntire:
			listeG4InGene = listeG4InGeneEntire[ l['idGene'] ]
			for G4InGene in listeG4InGene:
				startG4 = int(G4InGene.split('|')[1])
				endG4 = int(G4InGene.split('|')[2])
				coordG4 = [startG4,endG4]
				if feature == 'Transcript':
					g4inTranscript = isG4InTranscript(l['Strand'],
									coordG4, borderTranscript)
				elif feature == 'Junction':
					g4inTranscript = isJuncionInTranscript(coordG4,
											l['intronList'])
				annotationTranscript = AnnotationTranscript[ l['idTr'] ]
				if g4inTranscript and annotationTranscript: # PBM ICI COORD
					informationsOfG4 = list(G4DetectedInGene[G4InGene])
					if feature == 'Transcript':
						locationInTranscript = \
							getlocationG4InTranscript(BiotypeByTranscript,
							ProteinCoding, l['idTr'],
							coordG4, l['end5'], l['start3'],
							l['exonList'], l['intronList'],
							l['Strand'])
					elif feature == 'Junction':
						locationInTranscript = \
							getlocationG4InJunction(BiotypeByTranscript,
							ProteinCoding, l['idTr'])
					descriptionG4 = l['Chromosome'] + ':' + str(startG4) + '-' \
									+ str(endG4) + '|' + l['Strand']
					G4InTranscript = addG4InTranscriptome(G4InTranscript,
									l['idTr'],
									descriptionG4,
									informationsOfG4,
									locationInTranscript,
									biotypeTranscript)
					G4InGenome = addG4InGenome(G4InGenome,
												l['idGene'],
												descriptionG4,
												locationInTranscript)
					TranscriptPerG4 = addTranscriptPerG4(TranscriptPerG4,
									descriptionG4,l['idTr'])
	inputfile.close()
	return G4InTranscript, G4InGenome, TranscriptPerG4

def getlisteG4InGene(G4Detected):
	"""Searches all pG4 of a gene.

	:param G4Detected: all pG4 with its score and sequence.
	:type G4Detected: dictionary

	:returns: listeG4InGene, for a gene all pG4 in it.
	:rtype: list
	"""
	listeG4InGene = {}
	for descriptionG4 in G4Detected:
		geneID = descriptionG4.split('|')[0]
		if geneID not in listeG4InGene:
			listeG4 = []
		else:
			listeG4 = listeG4InGene.get(geneID)
		listeG4.append(descriptionG4)
		listeG4InGene[geneID] = listeG4
	return listeG4InGene

def returnG4InGene(G4DetectedInGene, inputfile, dicoParam, StrandByGene):
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
	:param StrandByGene: strand for each genes.
	:type StrandByGene: dictionary

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
			dicoLine = readLineG4Screener(line, StrandByGene, 'Gene')
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
						startG4 = dicoLine['GeneEnd'] - \
								(dicoLine['WindowStart'] -\
								dicoLine['GeneStart'])
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
					headerG4 = rF.createIdG4(dicoLine['GeneID'],
								startG4, endG4, dicoLine['Strand'])
					if headerG4 not in G4DetectedInGene and dicoLine['Strand']:
						G4DetectedInGene[headerG4] = [str(meanCGcC),
													str(meanG4Hunter),
													sequenceG4,
													str(meanG4NN)]
	inputfile.close()
	return G4DetectedInGene

def returnG4InJunction(G4DetectedInJunction, inputfile, dicoParam,
						StrandByGene):
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
	:param StrandByGene: strand for each genes.
	:type StrandByGene: dictionary

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
			dicoLine = readLineG4Screener(line, StrandByGene, 'Junction')
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
					startFirstWindow = dicoLine['JunctionStart']
					endFirstWindow = dicoLine['JunctionEnd']
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
					startG4 = getChromosomalPositionForJunction(startG4,
									dicoLine['Strand'], dicoParam['Window'],
									startFirstWindow, endFirstWindow)
					endG4 = getChromosomalPositionForJunction(endG4,
									dicoLine['Strand'], dicoParam['Window'],
									startFirstWindow, endFirstWindow)
					headerG4 = rF.createIdG4(dicoLine['GeneID'],
								startG4, endG4, dicoLine['Strand'])
					onJonction = isG4OnJunction(startG4, endG4,
								dicoLine['JunctionStart'],
								dicoLine['JunctionEnd'])
					if headerG4 not in G4DetectedInJunction and onJonction:
						G4DetectedInJunction[headerG4] = [str(meanCGcC),
														str(meanG4Hunter),
														sequenceG4,
														str(meanG4NN)]
	inputfile.close()
	return G4DetectedInJunction

def createDictionaryStrandByGene(filename):
	"""Retrieves strand for all genes of the chromosome.

	:param filename: path of the index file.
	:type filename: string

	:returns: dico, contains strand for all genes of the chromosome.
	:rtype: dictionary
	"""
	dico = {}
	inputfile = open(filename,'r')
	for line in inputfile:
		words = line.split('|')
		gene = words[1]
		strand = words[3]
		if gene not in dico:
			dico[gene]=strand
	inputfile.close()
	return dico

def createListCodingProtein():
	"""Creates a list all biotype from the coding protein.

	:returns: codingProtein, contain all biotype that are coding protein.
	:rtype: list
	"""
	codingProtein = ['IG_C_gene', 'IG_D_gene', 'IG_J_gene', 'IG_LV_gene',
					'IG_M_gene', 'IG_V_gene', 'IG_Z_gene',
					'nonsense_mediated_decay', 'protein_coding',
					'nontranslating_CDS', 'non_stop_decay',
					'TR_C_gene', 'TR_D_gene', 'TR_gene',
					'TR_J_gene', 'TR_V_gene']
	return codingProtein

def build_arg_parser():
	""" Retrieves paramters given in the command line.
	"""
	parser = argparse.ArgumentParser(description = 'G4Annotation')
	GITDIR=os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-chr', '--chromosome', default = 'X')
	parser.add_argument ('-specie', '--specie', default = 'HS')
	parser.add_argument ('-G4H', '--THRESHOLD_G4H', default = 0.9)
	parser.add_argument ('-CGCC', '--THRESHOLD_CGCC', default = 4.5)
	parser.add_argument ('-G4NN', '--THRESHOLD_G4NN', default = 0.5)
	parser.add_argument ('-E', '--EXTENSION', default = 100)
	parser.add_argument ('-W', '--WINDOW', default = 60)
	parser.add_argument ('-S', '--STEP', default = 10)
	return parser

def main(path, chromosome, specie, dicoParam):
	"""Main function.

	The main function contains 4 parts :

	- initialization of all lists and dictionnary containing informations
	- merge of all windows over thresholds. If two windows over thresolds
		are separated by an other window (under thresholds) between them, it
		be two differents pG4
	- annotation of pG4 : on genome and transcriptom
	- writing files

	"""
	GITDIR = os.getcwd()+'/'
	ProteinCoding = createListCodingProtein()
	G4DetectedInGene = {}
	G4DetectedInJunction = {}
	G4InTranscript = {}
	G4InGenome = {}
	TranscriptPerG4 = {}
	directory = path + '/chr' + chromosome
	# variable directory which contain the data for this chromosome
	index = directory + '/' + specie + '_transcript_unspliced_chr' + \
			chromosome + '_Index.txt'
	# file which contain info by transcript for this chromosome
	indexBiotypeTranscript = path + '/transcriptType/transcriptType_chr' + \
							chromosome
	# file which contain biotype of transcript for this chromosome
	print "Chromosome " + chromosome
	BiotypeByTranscript = \
		rF.createDictionaryBiotypeByTranscript(indexBiotypeTranscript)
	StrandByGene = createDictionaryStrandByGene(index)
	AnnotationTranscript = rF.GetAnnotationTranscript(index,
							ProteinCoding, BiotypeByTranscript)
	# get g4 from the ouput of G4RNA Screener
	for path, dirs, files in os.walk(directory):
		# for each element of the directory to passed
		for filename in files: # for each files (all csv)
			inputfile = directory + '/' + filename
			if ('gene_unspliced' in filename and '.csv' in filename ):
				G4DetectedInGene = returnG4InGene(G4DetectedInGene,
									inputfile, dicoParam,
									StrandByGene)
			elif ('Junction' in filename and '.fas' in filename):
				G4DetectedInJunction = returnG4InJunction(G4DetectedInJunction,
										inputfile, dicoParam,
										StrandByGene)
	listeG4InGeneEntire = getlisteG4InGene(G4DetectedInGene)
	listeG4InGeneJunction = getlisteG4InGene(G4DetectedInJunction)
	G4InTranscript, G4InGenome, TranscriptPerG4 = getInfoAboutpG4(index,
					BiotypeByTranscript, listeG4InGeneEntire,
					G4InTranscript, G4InGenome, TranscriptPerG4,
					AnnotationTranscript, G4DetectedInGene, ProteinCoding,
					'Transcript')
	G4InTranscript, G4InGenome, TranscriptPerG4 = getInfoAboutpG4(index,
					BiotypeByTranscript, listeG4InGeneJunction,
					G4InTranscript, G4InGenome, TranscriptPerG4,
					AnnotationTranscript, G4DetectedInJunction, ProteinCoding,
					'Junction')
	extractionG4InTranscript(GITDIR + 'results/perChromosome', specie,
							chromosome, G4InTranscript)
	extractionG4InGenome(GITDIR + 'results/perChromosome', specie, chromosome,
						G4InGenome)
	extractionTranscriptPerG4(GITDIR + 'results/perChromosome', specie,
							chromosome, TranscriptPerG4)
	print "\t Done."


if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path # directory containing all chromosome
	chromosome = arg.chromosome # chromosome to analyze
	specie = arg.specie # specie to analyse
	dicoParam = rF.createDicoParam(arg)
	main(path, chromosome, specie, dicoParam)
