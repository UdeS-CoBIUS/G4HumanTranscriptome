#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	March 2019

Description:
	This file is more like a librairy than a script. It contains all function
		that are used in 2 or more scripts. This allow to avoid redundance in
		other code.

"""

import pandas as pd
import re

def createDicoFamily():
	"""Creates a dictionnary with all classes and subclasses of transcripts.

	:returns: dicoFam, contains all classes and subclasses.
	:rtype: dictionary
	"""
	dicoFam = {'Coding' : ['IG_C_gene', 'IG_D_gene', 'IG_J_gene',
							'IG_LV_gene', 'IG_M_gene', 'IG_V_gene',
							'IG_Z_gene', 'nonsense_mediated_decay',
							'nontranslating_CDS', 'non_stop_decay',
							'protein_coding', 'TR_C_gene', 'TR_D_gene',
							'TR_gene', 'TR_J_gene', 'TR_V_gene'],
				'Pseudogene' : ['transcribed_unitary_pseudogene',
								'disrupted_domain', 'IG_C_pseudogene',
								'IG_J_pseudogene', 'IG_pseudogene',
								'IG_V_pseudogene', 'processed_pseudogene',
								'pseudogene',
								'transcribed_processed_pseudogene',
								'transcribed_unprocessed_pseudogene',
								'translated_processed_pseudogene',
								'translated_unprocessed_pseudogene',
								'TR_J_pseudogene', 'TR_V_pseudogene',
								'unitary_pseudogene', 'unprocessed_pseudogene',
								'polymorphic_pseudogene'],
				'LongNC' : ['macro_lncRNA', 'bidirectional_promoter_lncRNA',
							'sense_intronic', '3prime_overlapping_ncRNA',
							'ambiguous_orf', 'antisense',
							'lincRNA', 'ncrna_host','non_coding',
							'processed_transcript', 'retained_intron',
							'sense_overlapping'],
				'ShortNC' : ['vaultRNA', 'scaRNA', 'miRNA',
							'miRNA_pseudogene', 'misc_RNA',
							'misc_RNA_pseudogene', 'Mt_rRNA',
							'Mt_tRNA', 'Mt_tRNA_pseudogene',
							'ncRNA', 'pre_miRNA', 'RNase_MRP_RNA',
							'RNase_P_RNA', 'rRNA', 'rRNA_pseudogene',
							'scRNA', 'scRNA_pseudogene', 'snlRNA',
							'snoRNA', 'snoRNA_pseudogene', 'snRNA',
							'snRNA_pseudogene', 'SRP_RNA', 'tmRNA',
							'tRNA', 'tRNA_pseudogene','ribozyme'],
				'Predictif' : ['TEC']}
	return dicoFam

def createDicoFamilyFiltered():
	"""Creates a dictionnary with all classes and subclasses of transcripts.

	:returns: dicoFam, contains all classes and subclasses.
	:rtype: dictionary
	"""
	dicoFam = {'Coding' : ['nonsense_mediated_decay', 'non_stop_decay',
							'protein_coding'],
				'Pseudogene' : ['processed_pseudogene',
								'transcribed_unitary_pseudogene',
								'transcribed_processed_pseudogene',
								'transcribed_unprocessed_pseudogene',
								'translated_processed_pseudogene',
								'translated_unprocessed_pseudogene',
								'unitary_pseudogene', 'unprocessed_pseudogene'],
				'LongNC' : ['sense_intronic', 'antisense', 'lincRNA',
							'processed_transcript', 'retained_intron',
							'sense_overlapping'],
				'ShortNC' : ['vaultRNA', 'scaRNA', 'miRNA', 'misc_RNA', 'ncRNA',
							'pre_miRNA', 'scRNA', 'snlRNA', 'snoRNA' 'snRNA',
							'ribozyme'],
				'Predictif' : ['TEC']}
	return dicoFam

def createDictionaryBiotypeByTranscript(filename):
	"""Creates dictionary with transcripts subclasses.

	:param filename: name of file which contains the liste of biotype for each
		transcript of this chromosome.
	:type filename: string

	:returns: dico, contains transcripts and their subclasses.
	:rtype: dictionary
	"""
	dico = {}
	inputfile = open(filename,'r')
	for line in inputfile:
		if line != 'Gene stable ID\tTranscript stable ID\tChromosome/scaffold name\tTranscript type\n':
			words = line.split('\t')
			idTr = words[1]
			biotypeTranscript = words[3].rstrip()
			if idTr not in dico:
				dico[idTr] = biotypeTranscript
	inputfile.close()
	return dico

def GetAnnotationTranscript(filename, Coding, BiotypeByTranscript):
	"""Returns a dictionnary with the annotation quality of transcripts.

	In Ensembl, there is some transcripts with bad annotation. The ones we can
	checked are if a non coding transcripts have an UTR, which chould not exist.

	:param filename: name of the file which contains information for each
		transcript.
	:type filename: string
	:param Coding: contains all subclasses of coding transcripts.
	:type Coding: list
	:param BiotypeByTranscript: contains all transcripts and their biotypes.
	:type BiotypeByTranscript: dictionary

	:returns: AnnotationTranscript, contains all transcripts and if their
		annotation is good.
	:rtype: dictionary
	"""
	AnnotationTranscript = {}
	inputfile = open(filename,"r")
	for line in inputfile:
		words = line.split('|')
		transcriptId = words[0]
		start5 = words[7]
		end5 = words[8]
		start3 = words[9]
		end3 = words[10].rstrip()
		answer = True # variable answer by defaul True
		transcriptBiotype = ''
		transcriptBiotype = BiotypeByTranscript.get(transcriptId)
		if transcriptBiotype not in Coding:
			if (start5!='' or end5!='' or start3!='' or end3!='' ):
				# but if transcript has a 5'UTR or an 3' UTR bad anotation
				answer = False
		AnnotationTranscript[transcriptId] = answer
	inputfile.close()
	return AnnotationTranscript

def GetLengthFraction(positionA, positionB):
	"""Computes the length of a region.

	This fonction defines the lenght between two positions
	(positionA and positionB). It doesn't matter the order of position
	(positionA can be > positionB).

	:param positionA: first coordinate.
	:type positionA: string
	:param positionB: scond coordinates.
	:type positionB: string

	:returns: length, length of the region delimited by positionA and positionB.
	:rtype: integer
	"""
	length = 0
	if ((positionA and positionB) != ''):
		length = ( abs( int(positionA) - int(positionB) ) ) +1
	return length

def createDicoParam(arg):
	"""Creates a dictionary containing all parameters.

	:param arg: contains all parameters from G4RNA screener and the length
		of junctions.
	:type arg: arg

	:returns: dicoParam, contains all parameters from G4RNA screener
		and the length of junctions
	:rtype: dictionary
	"""
	dicoParam = {'g4H' : float(arg.THRESHOLD_G4H),
				'cGcC' : float(arg.THRESHOLD_CGCC),
				'g4NN' : float(arg.THRESHOLD_G4NN),
				'Extension' : int(arg.EXTENSION),
				'Window' : int(arg.WINDOW),
				'Step' : int(arg.STEP)}
	return dicoParam

def addWindowToG4Seq(g4Seq, windowSeq, step, windowLength):
	"""Adds a wwindow's sequence to the sequence of a pG4.

	:param g4Seq: sequence of a pG4.
	:type g4Seq: string
	:param windowSeq: sequence of the window.
	:type windowSeq: string
	:param step: step between two windows (parameter of G4RNA screener).
	:type step: integer
	:param windowLength: length of windows (parameter of G4RNA screener).
	:type windowLength: integer

	:returns: g4Seq, merge between the window and the old pG4.
	:rtype: string
	"""
	if len(windowSeq) < windowLength:
		# if the sequence is shorter than the sequence of the windows
		# (generraly last windows above the thresolds)
		g4Seq += windowSeq[-(len(windowSeq)-(windowLength-step)):]
	else: #
		g4Seq += windowSeq[-step:] # take the stepsize added
	return g4Seq

def iniDicoNucl():
	nuclCount = {'A':0,
				'T' : 0,
				'G' : 0,
				'C' : 0,
				'N' : 0}
	return nuclCount

def createIdG4(gene, startG4, endG4, strand):
	"""Creates the id of a pG4.

	:param gene: id of a gene.
	:type gene: string
	:param startG4: start of the pG4.
	:type startG4: integer
	:param endG4: end of the pG4.
	:type enG4: integer
	:param strand: strand of the pG4 and of the gene.
	:type strand: string

	:returns: the id of the pG4
	:rtype: string
	"""
	if strand == '1' or strand == '+':
		return gene + '|' + str(startG4) + '|' + str(endG4) + '|' + strand
	elif strand == '-1' or strand == '-':
		return gene + '|' + str(endG4) + '|' + str(startG4) + '|' + strand

def mean(liste):
   	return sum(liste)/len(liste)

def createIdpG4rShuffle(gene, startG4, endG4, strand):
	"""Creates the id of a pG4.

	:param gene: id of a gene.
	:type gene: string
	:param startG4: start of the pG4.
	:type startG4: integer
	:param endG4: end of the pG4.
	:type enG4: integer
	:param strand: strand of the pG4 and of the gene.
	:type strand: string

	:returns: the id of the pG4
	:rtype: string
	"""
	if strand == '1' or strand == '+':
		return gene + ';' + str(startG4) + ';' + str(endG4) + ';' + strand
	elif strand == '-1' or strand == '-':
		return gene + ';' + str(endG4) + ';' + str(startG4) + ';' + strand

def changeStrandFormat(strand):
	"""Changes the format of the strand from +/- to 1/-1.
	"""
	if strand == '+':
		strand = '1'
	elif strand == '-':
		strand = '-1'
	return strand

def addTypeTr(biotype):
	"""Apply function to retrieve the type of a transcript depending on biotype.

	:param biotype: biotype of a transcript.
	:type biotype: string

	:returns: type, coding or non coding.
	:rtype: string
	"""
	coding = ['IG_C_gene', 'IG_D_gene', 'IG_J_gene',
			'IG_LV_gene', 'IG_M_gene', 'IG_V_gene',
			'IG_Z_gene', 'nonsense_mediated_decay', 'non_stop_decay',
			'protein_coding', 'TR_C_gene', 'TR_D_gene',
			'TR_gene', 'TR_J_gene', 'TR_V_gene']
	if biotype in coding:
		return 'Coding'
	else:
		return 'Non coding'
