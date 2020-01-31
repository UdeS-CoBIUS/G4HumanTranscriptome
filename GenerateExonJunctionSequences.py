#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Contact:
	sarah.belhamiti@usherbrooke.ca

Date:
	December 2017

Description:
	This script aims to compute coordinates of introns and to create junctions
	sequences. Junctions corrsponds to the 100 upstreams nucleotides of an
	intron and to the 100 dowstream nucleotides of the same introns.
"""

import sys
import Bio
import os
import re
from pprint import pprint
from operator import itemgetter
import argparse

def importFastaChromosome(directory, chr):
    """Import the fasta file containing an entire chromosome sequence.

    :param directory: directory name containing all chromosome fasta.
    :type directory: string
    :param chr: chromosome wanted.
    :type chr: string

    :returns: sequence, entire chromosome sequence.
    :rtype: string
    """
    filename = directory + '/Homo_sapiens.GRCh38.dna.chromosome.' + \
                chr.split('r')[1] + '.fa'
    with open(filename) as f: # file opening
        content = f.read()
        l = content.split('\n')
        if l[0].startswith('>'):
            header = l[0]
            sequence = "".join(l[1:])
    return sequence

def readLineTranscriptInfo(l):
	"""Reads a line from the file containing infromations about transcripts.

	This file is downloaded from Ensembl via biomart. The human assembly is
	GRCh38p12.

	:param l: line from the transcript file.
	:type l: string

	:returns: dicoL, contains all info of the line but parsed.
	:rtype: dictionary
	"""
	words = l.split('\t')
	dicoL = {'geneId' : words[0],
			'trId' : words[1],
			'Chromosome' : words[2],
			'Biotype' : words[3],
			'5utrStart' : words[4],
			'5utrEnd' : words[5],
			'3utrStart' : words[6],
			'3utrEnd' : words[7],
			'exonStart' : words[8],
			'exonEnd' : words[9],
			'Rank' : words[10],
			'Strand' : words[11].rstrip()}
	return dicoL

def OrderInformationBiomart(directory, inputfilename):
	"""Creates one line per transcript with all feature in it.

	Here we read the file containing informations about transcripts and parse
	it. In the input file have informations about exons, UTR, chromosomes,
	strand and biotype. One line contains informations about one enxon. UTR are
	given on the first exon (5UTR) and on the last exon (3UTR). Introns coord
	are computated and added to the others. Then all informations are regrouped
	to fit in one line.

	:param directory: name of the directorycontaining the file.
	:type directory: string
	:param inputfilename: name of the file containing transcripts informations.
	:type inputfilename: string

	:returns: InformationPerGeneAndTranscript, corresponds to one line of the
		output.
	:rtype: string
	"""
	old_transcript = ""
	InformationPerGeneAndTranscript = {}
	inputfile = open(directory+inputfilename,"r")
	for line in inputfile:
		if 'ENS' in line:
			dicoL = readLineTranscriptInfo(line)
			if old_transcript != dicoL['trId']:
				var_rank = ""
				var_startExon = ""
				var_endExon = ""
				start5 = ""
				end5 = ""
				start3 = ""
				end3 = ""
				# if new transcript assignation
				var_rank = dicoL['Rank']
				var_startExon = dicoL['exonStart']
				var_endExon = dicoL['exonEnd']
				start5 = dicoL['5utrStart']
				end5 = dicoL['5utrEnd']
				start3 = dicoL['3utrStart']
				end3 = dicoL['3utrEnd']
				old_transcript = dicoL['trId']
			else: # add information for transcript with several exons
				# only 3UTR is update because it's given at the end of the tr
				var_rank = var_rank + ";" + dicoL['Rank']
				var_startExon = var_startExon + ";" + dicoL['exonStart']
				var_endExon = var_endExon + ";" + dicoL['exonEnd']
				start3 = dicoL['3utrStart']
				end3 = dicoL['3utrEnd']
			InformationPerGeneAndTranscript[dicoL['geneId']+"-"+dicoL['trId']] \
				= dicoL['Chromosome'] + "|" + dicoL['Biotype'] + "|" + \
				start5 + "|" + end5 + "|" + start3 + "|" + end3 + "|" + \
				var_startExon + "|" + var_endExon + "|" + var_rank + "|" + \
				dicoL['Strand']
	return InformationPerGeneAndTranscript

def ExonTotal(rank, strand, start_exon, end_exon):
	"""Retrieves from one transcript all its exons.

	:param rank: number of exon.
	:type rank: integer
	:param strand: strand of the transcript, can be 1 for forward and -1 for
		reverse.
	:type strand: string
	:param start_exon: all exon start for a transcript.
	:type start_exon: integer
	:param end_exon: all exon end for a transcript.
	:type end_exon: integer

	:returns: exon_total, array of array with start and end of each exon.
	:rtype: list
	"""
	i = 0
	exon_total = []
	while i < rank:
			exon = []
			if (strand == "1"): # if forward gene
				exon.append(int(start_exon[i]))
				exon.append(int(end_exon[i]))
				i += 1
			else: # if reverse gene
				exon.append(int(end_exon[i]))
				exon.append(int(start_exon[i]))
				i += 1
			exon_total.append(exon)
	return exon_total

def CreateStartIntron(exon_total, strand):
	"""Computes all introns starts for a transcript.

	:param exon_total: all start and end of each exon from a transcript.
	:type exon_total: list
	:param strand: strand of the transcript, can be 1 for forward and -1 for
		reverse.
	:type strand: string

	:returns: start_intron, all introns starts from a transcript.
	:rtype: list
	"""
	start_intron = []
	for i in exon_total:
		if i[1] != exon_total[-1][1]:
			# EXCEPT end last exon
			if strand == "1":
				intronStart = int( i[1] ) + 1
			else:
				intronStart = int( i[1] ) - 1
			start_intron.append(intronStart)
	return start_intron

def CreateEndIntron(exon_total, rank, strand):
	"""Computes all introns ends for a transcript.

	:param exon_total: all start and end of each exon from a transcript.
	:type exon_total: list
	:param rank: contain all exon ranks from a transcript.
	:type rank: list
	:param strand: strand of the transcript, can be 1 for forward and -1 for
		reverse.
	:type strand: string

	:returns: end_intron, all introns ends from a transcript.
	:rtype: list
	"""
	end_intron = []
	for i in exon_total:
		if int(i[0]) != int(exon_total[-len(rank)][0]):
			# EXCEPT start firt exon
			if strand == "1":
				introndEnd = int( i[0] ) - 1
			else:
				introndEnd = int( i[0] ) + 1
			end_intron.append(introndEnd)
	return end_intron

def IntronTotal(start_intron, end_intron):
	"""Gets all introns coordinates from a transcript.

	:param start_intron: all introns starts from a transcript.
	:type start_intron: list
	:param end_intron: all introns ends from a transcript.
	:type end_intron: list

	:returns: intron_total, all start and end of each intron from a transcript.
	:rtype: list
	"""
	i = 0
	intron_total = []
	while i < len(start_intron): # for number of intron
			intron = []
			intron.append( str(start_intron[i]) )
			intron.append( str(end_intron[i]) )
			intron_total.append(intron)
			i += 1
	return intron_total

def AddIntron(intron_total, Intron):
	""" Adds new introns in genes intron.

	:param intron_total: all introns from a transcript.
	:type intron_total: list
	:param Intron: all intron from transcriptome.
	:type Intron: list

	:returns: Intron, updated with new introns.
	:rtype: list
	"""
	for i in intron_total:
		if i not in Intron:
			Intron.append(i)
	return Intron

def AddTranscriptPerIntron(intron_total, Dico, transcriptID):
	"""Creates a dictionary with all introns and their transcript ID.

	:param intron_total: all introns from a transcript.
	:type intron_total: list
	:param Dico: contiains for each intron, the transcript it's belong.
	:type Dico: dictionary
	:param transcriptID: Ensembl transcript ID.
	:type transcriptID: string

	:returns: Dico, with new introns.
	:rtype: dictionary
	"""
	for i in intron_total:
		if i[0]+'-'+i[1] not in Dico:
			# if this intron is not in dico
			Dico[i[0]+'-'+i[1]] = transcriptID
		else:
			last_transcript = Dico.get( i[0]+'-'+i[1] )
			Dico[ i[0]+'-'+i[1] ] = last_transcript +"-"+ transcriptID
	return Dico

def reverseSequence(Sequence):
    """ Reverse complement a DNA sequence.

    :param Sequence: DNA sequence that will be reversed.
    :type Sequence: string

    :returns: Sequence, the initial DNA sequence but reverse complemented.
    :rtype: string
    """
    reverse = ""
    nucleotides = {'A' : 'T',
                    'T' : 'A',
                    'C' : 'G',
                    'G' : 'C'}
    for n in Sequence:
        if n in nucleotides:
            tmp = nucleotides[n]
        else :
            tmp = n # in some sequences there is many N or other letter
        reverse += tmp
    return reverse

def Fasta(chrSequence, start, end, strand):
	"""Gets a gene fasta sequence from Ensembl.

	:param geneID: Ensembl Id of a gene.
	:type geneID: string
	:param chromosome: name of the gene chromosome
	:type chromosome: string
	:param start: start of the gene.
	:type start: integer
	:param end: end of the gene.
	:type end: integer
	:param strand: strand of the gene.
	:type strand: string
	:param specie: name of the specie.
	:type specie: string

	:returns: DNA sequence of the gene.
	:rtype: string
	"""
	seq = chrSequence[start-1:end-1]
	# if strand == '1':
	#
	# else:
	# 	seq = chrSequence[end-1:start-1]
	# 	seq = reverseSequence(seq)
	return seq

def CreateSequence(dirFasta, dirData, inputfilename, IntronPerGene, InfoPerGene,
					extension, chromosome):
	"""Writes output with all junction sequences.

	:param directory: name of the directory containing the file.
	:type directory: string
	:param inputfilename: name of the file containing transcripts informations.
	:type inputfilename: string
	:param IntronPerGene: contains all genes and all its introns.
	:type IntronPerGene: dictionary
	:param InfoPerGene: contains informations of genes.
	:type InfoPerGene: dictionary
	:param extension: length of the upstream sequence and dowstream sequence we
		retrieve.
	:type extension: integer
	:param chromosome: current chromosome.
	:type chromosome: string
	"""
	chrSequence = importFastaChromosome(dirFasta, chromosome)
	output = open(dirData+inputfilename.split(".")[0]+"_Sequence.txt","w")
	for gene in IntronPerGene:# for every intron by gene
		if gene in InfoPerGene:
			geneID = gene
			chromosome = InfoPerGene[gene].split("\t")[0]
			biotype = InfoPerGene[gene].split("\t")[1]
			strand = InfoPerGene[gene].split("\t")[2]
			for intron in IntronPerGene[gene]:
				start_intron = int(intron[0])
				end_intron = int(intron[1])
				if strand == "1":
					start_sup = start_intron - 1 - extension
					end_sup = start_intron - 1
					start_inf = end_intron + 1
					end_inf = end_intron + 1 + extension
					sequence_amont = Fasta(chrSequence, start_sup, end_sup, strand)
					sequence_aval = Fasta(chrSequence, start_sup, end_sup,	strand)
				else :
					start_sup = start_intron + 1 + extension
					end_sup = start_intron + 1
					start_inf = end_intron - 1
					end_inf = end_intron - 1 - extension
					sequence_amont = Fasta(chrSequence, end_sup, start_sup, strand)
					sequence_aval = Fasta(chrSequence, end_sup, start_sup,	strand)
					sequence_amont = reverseSequence(sequence_amont)
					sequence_aval = reverseSequence(sequence_aval)
				sequence = sequence_amont + sequence_aval
				output.write(">" + geneID + "|" + "-".join(intron) + \
					"\n" + sequence + "\n")
	output.close()

def CreateIndex(directory, inputfilename, InfoPerGeneAndTranscript,
	ExonPerTranscript, IntronPerTranscript):
	"""Writes output with all transcripts informations.

	:param directory: name of the directorycontaining the file.
	:type directory: string
	:param inputfilename: name of the file containing transcripts informations.
	:type inputfilename: string
	:param InfoPerGeneAndTranscript: contains all genes and transcripts info.
	:type InfoPerGeneAndTranscript: dictionary
	:param ExonPerTranscript: contains all transcripts and their exons.
	:type ExonPerTranscript: dictionary
	:param IntronPerTranscript: contains all transcripts and all its introns.
	:type IntronPerTranscript: dictionary
	"""
	output = open(directory + inputfilename.split(".")[0] + "_Index.txt","w")
	for gene in InfoPerGeneAndTranscript:
		exonList = ""
		intronList = ""
		geneID = gene.split("-")[0]
		transcriptID = gene.split("-")[1]
		chromosome = InfoPerGeneAndTranscript[gene].split("|")[0]
		biotype = InfoPerGeneAndTranscript[gene].split("|")[1]
		start5 = InfoPerGeneAndTranscript[gene].split("|")[2]
		end5 = InfoPerGeneAndTranscript[gene].split("|")[3]
		start3 = InfoPerGeneAndTranscript[gene].split("|")[4]
		end3 = InfoPerGeneAndTranscript[gene].split("|")[5]
		start_exon = InfoPerGeneAndTranscript[gene].split("|")[6].split(";")
		strand = InfoPerGeneAndTranscript[gene].split("|")[9]
		if transcriptID in ExonPerTranscript:
				exon = ExonPerTranscript[transcriptID]
				for couple in exon:
					exonList = exonList+'-'.join(str(x) for x in couple)+";"
		else:
			print "error exon", geneID, transcriptID
		if transcriptID in IntronPerTranscript:
				intron = IntronPerTranscript[transcriptID]
				for couple in intron:
					intronList = intronList + '-'.join(couple) + ";"
		else:
			print "error intron", geneID, transcriptID
		output.write(transcriptID +"|"+ geneID +"|"+ chromosome +"|"+ strand +\
			"|"+ biotype +"|"+ exonList[:-1] +"|"+ intronList[:-1] +"|"+ \
			start5 +"|"+ end5 +"|"+ start3 +"|"+ end3 +"\n")
	output.close()

def build_arg_parser():
	GITDIR = os.getcwd()+'/'
	parser = argparse.ArgumentParser(description = 'ReplaceInformationBiomart')
	parser.add_argument ('-p', '--path', default = GITDIR)
	parser.add_argument ('-chr', '--chromosome', default = 'chrY')
	parser.add_argument ('-ext', '--extension', default = 100)
	return parser

def main(path, inputfilename, extension, chromosome):
	old_gene = ""
	IntronPerGene = {}
	ExonPerTranscript = {}
	IntronPerTranscript = {}
	TranscriptPerIntron = {}
	InfoPerGene = {}
	dirData = path + 'Data/' + chromosome + '/'
	dirFasta = path + 'Data/Fasta'
	InfoPerGeneAndTr = OrderInformationBiomart(dirData, inputfilename)
	for geneTr in InfoPerGeneAndTr:
		infoTr = InfoPerGeneAndTr[geneTr].split("|")
		geneID = geneTr.split("-")[0]
		transcriptID = geneTr.split("-")[1]
		chr = infoTr[0]
		biotype = infoTr[1]
		start5 = infoTr[2]
		end5 = infoTr[3]
		start3 = infoTr[4]
		end3 = infoTr[5]
		start_exon = infoTr[6].split(";")
		end_exon = infoTr[7].split(";")
		rank = infoTr[8].split(";")
		strand = infoTr[9]
		# if (re.compile('[1-9, X, Y]').search(chr) and not
		# 	re.compile('[A-W, Z]').search(chr)):
		InfoPerGene[geneID] = chr +'\t'+ biotype +'\t'+ strand
		exon_total = ExonTotal(len(rank), strand, start_exon, end_exon)
		# in some case, exon are not in the good order, so we sort them to
		# avoid to get bad introns (their coords ar computed using the order
		# of lines in the input file)
		if strand == '1':
			exon_total = sorted(exon_total, key=itemgetter(0))
		else:
			exon_total = sorted(exon_total, key=itemgetter(0), reverse=True)
		start_intron = CreateStartIntron(exon_total, strand)
		end_intron = CreateEndIntron(exon_total, rank, strand)
		intron_total = IntronTotal(start_intron, end_intron)
		ExonPerTranscript[transcriptID] = exon_total
		IntronPerTranscript[transcriptID] = intron_total
		if (old_gene != geneID):
			Intron = []
			Intron = AddIntron(intron_total, Intron)
			old_gene = geneID
		else:
			# if other transcript
			Intron = AddIntron(intron_total, Intron)
			old_gene = geneID
		IntronPerGene[geneID] = Intron
		TranscriptPerIntron = AddTranscriptPerIntron(intron_total,
			TranscriptPerIntron, transcriptID)
	##Create file Sequence
	CreateSequence(dirFasta, dirData, inputfilename, IntronPerGene, InfoPerGene,
		extension, chromosome)
	CreateIndex(dirData, inputfilename, InfoPerGeneAndTr, ExonPerTranscript,
		IntronPerTranscript)

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	chromosome = arg.chromosome
	extension = arg.extension
	inputfilename = 'HS_transcript_unspliced_'+chromosome+'.txt'
	main(path, inputfilename, extension, chromosome)
