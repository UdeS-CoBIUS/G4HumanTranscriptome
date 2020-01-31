#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
    Copyright Universite of Sherbrooke, departement of biochemistry and
    departement    of computation.

Date:
    June 2019

Description:
    This script aims to generate all random fasta needed. It will generate one
    fasta file for each segment location but also for junction. The id of each
    fasta is like : biotype:Trid:Start-End:Strand, exept for intron where the
    biotype is at the end. To generate random sequences we retrieve the sequence
    of the location we want, we transformed it into a list, then we shufle it
    (change all position of all element randomly), to finally join the shuffled
    list with no characters to recreate the sequence.

Data availability:
    * fasta sequences used with G4RNA screener were download from Ensembl. The human assembly is GRCh38p12.
    * transcripts information were downloaded from ensembl via ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/

Abbreviation:
    * bt for biotype
    * tr for transcript
    * dico for dictionnary
    * w for word or words
    * dir for directory
    * cpt for compteur (counter)
"""

import os
import re
import math
import random
import argparse
from Bio import SeqIO
from pprint import pprint
import recurentFunction as rF

def shuffleSeq(seq):
    """Shuffle a sequence.

    This function aims to shuffle a fasta sequence to randomize its sequence.
    The fasta sequence is imported from a fasta file, then converted into a
    list (one element corrresponds to one nucleotide), the list is shuffled and
    tehn joined with nothing to recreate the sequence.

    :param seq: sequence to shuffle.
    :type seq: string

    :returns: seq, sequence shuffled
    :rtype: string
    """
    seq = list(seq)
    random.shuffle(seq)
    seq = ''.join(seq)
    return seq

def importFastaChromosome(directory, chr):
    """Import the fasta file containing an entire chromosome sequence.

    :param directory: directory name containing all chromosome fasta.
    :type directory: string
    :param chr: chromosome wanted.
    :type chr: string

    :returns: sequence, entire chromosome sequence.
    :rtype: string
    """
    filename = directory + 'Homo_sapiens.GRCh38.dna.chromosome.' + str(chr) + '.fa'
    with open(filename) as f: # file opening
        content = f.read()
        l = content.split('\n')
        if l[0].startswith('>'):
            header = l[0]
            sequence = "".join(l[1:])
    return sequence

def writeFasta(outDir, dicoFasta, location):
    """From a dictionary {id : seq}, write a fasta file.

    :param outDir: name of the directory where the output file need to be writen.
    :type outDir: string
    :param dicoFasta: {id : seq}
    :type dicoFasta: dictionary
    :param location: name of the chromosome.
    :type location: string
    """
    output = open(outDir +'Shuffle_chr'+chr+'.fas', "w")
    for id in dicoFasta:
        output.write(id + "\n")
        nbLine = math.ceil( float( len(dicoFasta[id]) ) / 60 )
        cpt1 = 0
        cpt2 = 60
        for i in range(0,int(nbLine)) :
            output.write(dicoFasta[id][cpt1:cpt2] + "\n")
            # to have a new line after 60 characters
            cpt1 += 60
            cpt2 += 60
    output.close()

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

def getUTRCoordByStrand(w, dicoLine):
	"""Retrieves coords of UTR depending on the strand.

	:param w: contains all words of one line of the index file.
	:type w: list
	:param dicoLine: contains all info on a line (except UTR).
	:type dicoLine: dictionary

	:returns: dicoLine, contains all info on a line, updated with UTR.
	:rtype: dictionary
	"""
	if dicoLine['Strand'] == '1':
		dicoLine.update({'start5' : w[7],
						'end5' : w[8],
						'start3' : w[9],
						'end3' : w[10]})
	else:
		dicoLine.update({'start5' : w[8],
						'end5' : w[7],
						'start3' : w[10],
						'end3' : w[9]})
	return dicoLine

def parseIndex(w):
    """Transforme a line from the index file into a dictionary.

    As one line of the index file contain many information, we parsed it into
    a dictionnary to access easly to different informations.

    :param w: contains all words of one line of the index file.
	:type w: list

    :returns: dicoLine, contains all information of one index line.
    :rtype: dictionary
    """
    dicoLine = {'idTr' : w[0],
				'idGene' : w[1],
				'Chromosome' : w[2],
				'Strand' : w[3],
				'geneBiotype' : w[4],
				'exonList' : w[5].split(';'),
				'intronList' : w[6].split(';')}
    dicoLine.update(getUTRCoordByStrand(w, dicoLine))
    return dicoLine

def createCoordSpliceSite(coord1, coord2, dicoLine):
    """Depending on the strand, create coordinates of splice sites.

    Creates coordnates of Donor splice site and Acceptor splice site, depending
    of the strand. The coordinates are +-59 nucletides of the intron coords
    because we only want to detect G4 that are close to it.

    :param coord1: start of the intron.
	:type coord1: int
    :param coord2: end of the intron.
	:type coord2: int
    :param dicoLine: line of the index file.
	:type dicoLine: dictionary

    :returns: idDonor, idAcceptor, id of the Donor splice site and of the
        Acceptor splice sites.
    :rtype: strings
    """
    if dicoLine['Strand'] == '1':
        donorStart = coord1 - 59
        donorEnd = coord1 + 58
        acceptorStart = coord2 - 59
        acceptorEnd = coord2 + 58
    else:
        donorStart = coord1 + 59
        donorEnd = coord1 - 58
        acceptorStart = coord2 + 59
        acceptorEnd = coord2 - 58
    idDonor = 'donor:'+dicoLine['Chromosome']+':'+str(donorStart)+'-'+str(donorEnd)+':'+dicoLine['Strand']
    idAcceptor = 'acceptor:'+dicoLine['Chromosome']+':'+str(acceptorStart)+'-'+str(acceptorEnd)+':'+dicoLine['Strand']
    return idDonor, idAcceptor

def addExonIntron(type, dico, dicoLine):
    """Add to the index of location exon, intron and also spliced sites.

    First the exon or intron location is created and added to the index
    location. If the location is an intron, we create Donor spliced site and
    Acceptor splice site using the +- 59 nucleotides upstream and dowstream
    of each intron coords.

    :param type: name of the location, exon or intron.
	:type type: string
    :param dico: aims to contain all location of a chromosome.
	:type dico: dictionary
    :param dicoLine: line of the index file.
	:type dicoLine: dictionary

    :returns: dico, updated with new locations.
    :rtype: dictionary
    """
    listName = type+'List'
    for elem in dicoLine[listName]:
        if elem:
            elemId = type+':'+dicoLine['Chromosome']+':'+elem+':'+dicoLine['Strand']
            if elemId not in dico[ dicoLine['idGene'] ]:
                dico[ dicoLine['idGene'] ][elemId] = {}
            if type == 'intron' and elem:
                coord1 = int(elem.split('-')[0])
                coord2 = int(elem.split('-')[1])
                idDonor,idAcceptor = createCoordSpliceSite(coord1, coord2, dicoLine)
                if idDonor not in dico[ dicoLine['idGene'] ]:
                    dico[ dicoLine['idGene'] ][idDonor] = {}
                if idAcceptor not in dico[ dicoLine['idGene'] ]:
                    dico[ dicoLine['idGene'] ][idAcceptor] = {}
            if dicoLine['idTr'] in dicoBt:
                dico[ dicoLine['idGene'] ][elemId].update({dicoLine['idTr'] : dicoBt[ dicoLine['idTr'] ]})
                if type == 'intron' and elem:
                    dico[ dicoLine['idGene'] ][idDonor].update({dicoLine['idTr'] : dicoBt[ dicoLine['idTr'] ]})
                    dico[ dicoLine['idGene'] ][idAcceptor].update({dicoLine['idTr'] : dicoBt[ dicoLine['idTr'] ]})
    return dico

def overlaps(interval1, interval2):
    """Computes the distance or overlap between two interval.

    Compute : min(ends) - max(starts). If > 0, return the number of bp of
    overlap, if 0, they are book-ended and if < 0 return the distance in
    bp between them.

    :param interval1: contain the start and the end of a OQs.
    :type interval1: list of int
    :param interval2: contain the start and the end of a gff feature.
    :type interval2: list of int

    :returns: min(ends) - max(starts)
    :rtype: int
    """
    return min(interval1[1], interval2[1]) - max(interval1[0], interval2[0])

def findCDS(coord5UTR, coord3UTR, coordExon):
    """Finds where is the CDS in an exon.

    We compute the overlap between UTR and exon. Then, depending of the overlap,
    we find CDS coords. If the there is a perfect overlap UTR-Exon then there is
    no CDS. If there no overlap then the exon is an CDS. Finally if there is an
    overlap, we compute the CDS coords.

    :param coord5UTR: coords of a 5UTR.
	:type coord5UTR: list of int
    :param coord3UTR: coords of a 3UTR.
	:type coord3UTR: list of int
    :param coordExon: coords of an exon.
	:type coordExon: list of int

    :returns: coordCDS, coords of a CDS.
    :rtype: list of int
    """
    overlap5UTR = overlaps(coord5UTR, coordExon)
    overlap3UTR = overlaps(coord3UTR, coordExon)
    if (overlap5UTR == coord5UTR[1]-coord5UTR[0] or
    overlap3UTR == coord3UTR[1]-coord3UTR[0] or
    overlap5UTR == coord5UTR[0]-coord5UTR[1] or
    overlap3UTR == coord3UTR[0]-coord3UTR[1]):
        coordCDS = []
    elif overlap5UTR < 0 and overlap3UTR < 0:
        #Aucun overlap
        coordCDS = coordExon
    elif overlap5UTR > 0:
        if coordExon[0] < coordExon[1]:
            #forward strand
            startCDS = coord5UTR[1] + 1
            endCDS = coordExon[1]
            coordCDS = [startCDS, endCDS]
        else:
            #reverse strand
            startCDS = coord5UTR[1] - 1
            endCDS = coordExon[1]
            coordCDS = [startCDS, endCDS]
    elif overlap3UTR > 0:
        if coordExon[0] < coordExon[1]:
            #forward strand
            startCDS = coordExon[0]
            endCDS = coord3UTR[0] - 1
            coordCDS = [startCDS, endCDS]
        else:
            #reverse strand
            startCDS = coordExon[0]
            endCDS = coord3UTR[1] + 1
            coordCDS = [startCDS, endCDS]
    else:
        if coord5UTR[0] == coord5UTR[1] or coord3UTR[0] == coord3UTR[1]:
            coordCDS = [ coordExon[0], coordExon[1] ]
    return coordCDS

def addUTRCDS(dicoLine, dico, dicoBt):
    """Add to the index of location UTR and CDS.

    CDS coords are computed from the comparison of exon coords and UTR coords.
    All those information are added to the index of location.

    :param dicoLine: contains all info on a line (except UTR).
	:type dicoLine: dictionary
    :param dico: aims to contain all location of a chromosome.
	:type dico: dictionary
    :param dicoBt: contains all biotype of all transcripts.
	:type dicoBt: dictionary

    :returns: dico, updated with new locations.
    :rtype: dictionary
    """
    exonList = [ [int(e.split('-')[0]), int(e.split('-')[1])] for e in dicoLine['exonList'] ]
    if dicoLine['start5'] and dicoLine['end5']:
        coords5UTr = [ int(dicoLine['start5']), int(dicoLine['end5']) ]
        id5UTR = '5UTR:'+dicoLine['Chromosome']+':'+dicoLine['start5']+'-'+dicoLine['end5']+':'+dicoLine['Strand']
        if id5UTR not in dico[ dicoLine['idGene'] ]:
            dico[ dicoLine['idGene'] ][id5UTR] = {}
        if dicoLine['idTr'] in dicoBt:
            dico[ dicoLine['idGene'] ][id5UTR].update({dicoLine['idTr'] : dicoBt[ dicoLine['idTr'] ]})
    else:
        coords5UTr = [0, 0]
    if dicoLine['start3'] and dicoLine['end3']:
        coords3UTR = [ int(dicoLine['start3']), int(dicoLine['end3']) ]
        id3UTR = '3UTR:'+dicoLine['Chromosome']+':'+dicoLine['start3']+'-'+dicoLine['end3']+':'+dicoLine['Strand']
        if id3UTR not in dico[ dicoLine['idGene'] ]:
            dico[ dicoLine['idGene'] ][id3UTR] = {}
        if dicoLine['idTr'] in dicoBt:
            dico[ dicoLine['idGene'] ][id3UTR].update({dicoLine['idTr'] : dicoBt[ dicoLine['idTr'] ]})
    else:
        coords3UTR = [0, 0]
    if dicoLine['idTr'] in dicoBt:
        if rF.addTypeTr(dicoBt[ dicoLine['idTr'] ]) == 'Coding':
            coordsCDS = [ findCDS(coords5UTr, coords3UTR, e) for e in  exonList]
            coordsCDS = [CDS for CDS in coordsCDS if CDS]
            for CDS in coordsCDS:
                    CDS = '-'.join(str(x) for x in CDS)
                    elemId = 'CDS:'+dicoLine['Chromosome']+':'+CDS+':'+dicoLine['Strand']
                    if elemId not in dico[ dicoLine['idGene'] ]:
                        dico[ dicoLine['idGene'] ][elemId] = {}
                    if dicoLine['idTr'] in dicoBt:
                        dico[ dicoLine['idGene'] ][elemId].update({dicoLine['idTr'] : dicoBt[ dicoLine['idTr'] ]})
    return dico

def getDicoIndex(filename, dicoBt):
    """Parses the Index file into a dictionary to get all locations.

    Browse all transcript in the index file and computes CDS, Donor splice site
    and Acceptor splice site coordinates from other coords feature (exon,
    intron and UTR).

    :param filename: name of the index file,
	:type filename: string
    :param dicoBt: contains bt of all transcripts.
	:type dicoBt: dictionary

    :returns: dico, index file parsed to obtaine all locations informations.
    :rtype: dictionary
    """
    dico = {}
    with open(filename) as f:
        lines = f.read().splitlines()
        for l in lines:
            l = l.rstrip()
            words = l.split('|')
            dicoLine = parseIndex(words)
            if dicoLine['idGene'] not in dico:
                dico[ dicoLine['idGene'] ] = {}
            dico.update(addExonIntron('exon', dico, dicoLine))
            dico.update(addExonIntron('intron', dico, dicoLine))
            if dicoLine['start5'] or dicoLine['start3']:
                dico.update(addUTRCDS(dicoLine, dico, dicoBt))
    return dico

def createFasta(dico, fastaFile, chr, outputDir):
    """Create fasta file from the index of location.

    Each sequence corresponds to one location. The id is like :
    Gene:location:chr:start-end:Tr1-bt|Tr2-bt
    Location of only one nucleotides are not added to the fasta file.
    The sequence is retrieved from a fasta file containing the entire
    chromosome. This file have been dowloaded from the ensembl FTP.

    :param dico: aims to contain all location of a chromosome.
	:type dico: dictionary
    :param fastaFile: directory containing fasta of all chromosomes.
	:type fastaFile: string
    :param chr: name of the chromosome.
	:type chr: string
    :param outputDir: name of the output directory.
	:type outputDir: string
    """
    fastaRandom = {}
    chrSeq = importFastaChromosome(fastaFile, chr)
    for gene in dico:
        for location in dico[gene]:
            start = int(location.split(':')[2].split('-')[0])
            end = int(location.split(':')[2].split('-')[1])
            if end != start:
                strand = location.split(':')[3]
                if strand == '1':
                    seq = chrSeq[start-1:end-1]
                else:
                    seq = chrSeq[end-1:start-1]
                    seq = reverseSequence(seq)
                randomSeq = shuffleSeq(seq)
                listTr = [ tr+'-'+dico[gene][location][tr] for tr in dico[gene][location] ]
                id = '>'+gene+':'+location+':'+'|'.join(listTr)
                fastaRandom[id] = randomSeq
    writeFasta(outputDir, fastaRandom, chr)

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'generateRandom')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    parser.add_argument ('-chr', '--chromosome', default = 'Y')
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    chr = arg.chromosome
    path = arg.path
    fastaFile = path+'Fasta/'
    outputDir = path+'/chr'+chr+'/'
    print(chr)
    dicoBt = rF.createDictionaryBiotypeByTranscript(path+
            '/transcriptType/transcriptType_chr'+chr)
    dicoInfo = getDicoIndex(path+'/chr'+chr+'/HS_transcript_unspliced_chr'+
            chr+'_Index.txt', dicoBt)
    createFasta(dicoInfo, fastaFile, chr, outputDir)
