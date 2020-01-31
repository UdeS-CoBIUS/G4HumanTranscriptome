#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	This script filter pG4r and rG4 fasta files, depending on their comon
    transcript ID. pG4 found in non mRNA transcripts, or on location that do not
    exist on mRNA are removed.

Usage:
    python ~/PATH/FilterG4DependingOnTr.py
"""

import re
import sys
import os
import math
import argparse
from Bio import SeqIO
from pprint import pprint

def writeFasta(outFN, dicoFasta):
    """Writes a fasta file from a dictionary.

    :param outFN: FileName of the output.
	:type outFN: string
    :param dicoFasta: sequences id and sequences.
	:type dicoFasta: dictionary

    :returns: pG4Filtered
    :rtype: dictionary
    """
    output = open(outFN, "w")
    for id in dicoFasta:
        output.write('>' + id + "\n")
        nbLine = math.ceil( float( len(dicoFasta[id]) ) / 60 )
        cpt1 = 0
        cpt2 = 60
        for i in range(0,int(nbLine)) :
            output.write(dicoFasta[id][cpt1:cpt2] + "\n")
            # to have a new line after 60 characters
            cpt1 += 60
            cpt2 += 60
    output.close()

def importpG4seq(filename, IdDico, listTrrG4):
    """Imports pG4r dataset.

    This function aims to import pG4r in a dictionary. pG4r in introns are
    removed.

    :param filename: name of the file with rG4.
	:type filename: string
    :param IdDico: correspondance between ensembl and refseq ids.
	:type IdDico: dictionary
    :param listTrrG4: list of transcripts with rG4.
	:type listTrrG4: list

    :returns: pG4Filtered
    :rtype: dictionary
    """
    tmpDico = {}
    fastaOrigin = SeqIO.parse(open(filename),'fasta')
    pG4Filtered = {}
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        name = name.split('|')
        location = name[3]
        if name[0] in IdDico['Ensembl to Refseq coding'] and name[0] in listTrrG4 and 'Intron' not in location:
            pG4Filtered[ ':'.join(name) ] = seq
            if location not in tmpDico:
                tmpDico[location] = 1
            else:
                tmpDico[location] += 1
        elif name[0] in IdDico['Ensembl to Refseq non coding'] and name[0] in listTrrG4  and 'Intron' not in location:
            pG4Filtered[ ':'.join(name) ] = seq
            if location not in tmpDico:
                tmpDico[location] = 1
            else:
                tmpDico[location] += 1
    pprint(tmpDico)
    return pG4Filtered

def importRG4seqPDS(filename, IdDico):
    """Imports K+PDS rG4 dataset.

    This function aims to import rG4 in a dictionary. It also creates a list of
    transcript with rG4. rG4 are filtered depending on their sequences, if a
    minimal motif G(2+)X(+)G(2+)X(+)G(2+)X(+)G(2+) is not found then we don't
    keep the rG4. It also prints a dictionary with the number of rG4
    by locations.

    :param filename: name of the file with rG4.
	:type filename: string
    :param IdDico: correspondance between ensembl and refseq ids.
	:type IdDico: dictionary

    :returns: rG4Filtered, listTrrG4
    :rtype: dictionary and list
    """
    listTrrG4 = []
    tmpDico = {}
    minMotif = r"([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,})"
    fastaOrigin = SeqIO.parse(open(filename),'fasta')
    rG4Filtered = {}
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        name = name.split(':')
        listG = re.findall(minMotif, seq)
        location = name[5]
        if listG:
            if name[4] in IdDico['Refseq coding to Ensembl'] and 'Intron' not in location:
                name[4] = IdDico['Refseq coding to Ensembl'][ name[4] ]
                rG4Filtered[ ':'.join(name) ] = seq
                listTrrG4.append(name[4])
                if location not in tmpDico:
                    tmpDico[location] = 1
                else:
                    tmpDico[location] += 1
            elif name[4] in IdDico['Refseq non coding to Ensembl'] and 'Intron' not in location:
                name[4] = IdDico['Refseq non coding to Ensembl'][ name[4] ]
                rG4Filtered[ ':'.join(name) ] = seq
                listTrrG4.append(name[1])
                if location not in tmpDico:
                    tmpDico[location] = 1
                else:
                    tmpDico[location] += 1
    pprint(tmpDico)
    return rG4Filtered, listTrrG4

def importRG4seq(filename, IdDico):
    """Imports K+ rG4 dataset.

    This function aims to import rG4 in a dictionary. It also creates a list of
    transcript with rG4. rG4 are filtered depending on their sequences, if a
    minimal motif G(2+)X(+)G(2+)X(+)G(2+)X(+)G(2+) is not found then we don't
    keep the rG4. It also prints a dictionary with the number of rG4
    by locations.

    :param filename: name of the file with rG4.
	:type filename: string
    :param IdDico: correspondance between ensembl and refseq ids.
	:type IdDico: dictionary

    :returns: rG4Filtered, listTrrG4
    :rtype: dictionary and list
    """
    listTrrG4 = []
    tmpDico = {}
    minMotif = r"([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,}[A-Z,a-z]{1,})([G,g]{2,})"
    fastaOrigin = SeqIO.parse(open(filename),'fasta')
    rG4Filtered = {}
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        name = name.split(':')
        listG = re.findall(minMotif, seq)
        location = name[6]
        if listG:
           if name[1] in IdDico['Refseq coding to Ensembl'] and 'Intron' not in location:
               name[1] = IdDico['Refseq coding to Ensembl'][ name[1] ]
               rG4Filtered[ ':'.join(name) ] = seq
               listTrrG4.append(name[1])
               if location not in tmpDico:
                   tmpDico[location] = 1
               else:
                   tmpDico[location] += 1
           elif name[1] in IdDico['Refseq non coding to Ensembl'] and 'Intron' not in location:
               name[1] = IdDico['Refseq non coding to Ensembl'][ name[1] ]
               rG4Filtered[ ':'.join(name) ] = seq
               listTrrG4.append(name[1])
               if location not in tmpDico:
                   tmpDico[location] = 1
               else:
                   tmpDico[location] += 1
    pprint(tmpDico)
    return rG4Filtered, listTrrG4

def getIdTr(filename):
    """Gets a dictionary of Tr id correspondance.

    This function aims to create a dictionary that make correrspondance between
    Ensembl id and refseq id.

    :param filename: name of the file with transcript id correspondance.
	:type filename: string

    :returns: IdDico, contains correrspondance between Ensembl id and refseq id.
    :rtype: dictionary
    """
    IdDico = {'Ensembl to Refseq coding': {},
            'Ensembl to Refseq non coding': {},
            'Refseq coding to Ensembl': {},
            'Refseq non coding to Ensembl': {}}
    with open(filename) as f:
        content = f.read()
        lines = content.split('\n')
        del lines[0]
        del lines[-1]
        for l in lines:
            w = l.split('\t')
            if w[1]:
                if w[0] not in IdDico['Ensembl to Refseq coding']:
                    IdDico['Ensembl to Refseq coding'][w[0]] = w[1]
                if w[1] not in IdDico['Refseq coding to Ensembl']:
                    IdDico['Refseq coding to Ensembl'][w[1]] = w[0]
            if w[2]:
                if w[0] not in IdDico['Ensembl to Refseq non coding']:
                    IdDico['Ensembl to Refseq non coding'][w[0]] = w[2]
                if w[2] not in IdDico['Refseq non coding to Ensembl']:
                    IdDico['Refseq non coding to Ensembl'][w[2]] = w[0]
    return IdDico

def main(path):
    IdFN = path + 'Data/rG4seq/listTr/HS_transcript_ID.txt'
    rG4FastaFN = path + 'Data/rG4seq/cdt_K.fas'
    rG4PDSFastaFN = path + 'Data/rG4seq/cdt_PDS-K.fas'
    pG4FastaFN = path + 'Results/All/HS_All_G4InTranscript.fa'
    IdDico = getIdTr(IdFN)
    dicoFastarG4, listTrrG4 = importRG4seq(rG4FastaFN, IdDico)
    dicoFastarG4PDS, listTrrG4 = importRG4seqPDS(rG4PDSFastaFN, IdDico)
    dicoFastapG4 = importpG4seq(pG4FastaFN, IdDico, listTrrG4)
    writeFasta(path + 'Data/rG4seq/rG4_Filtered.fas', dicoFastarG4)
    writeFasta(path + 'Data/rG4seq/rG4PDS_Filtered.fas', dicoFastarG4PDS)
    writeFasta(path + 'Data/rG4seq/pG4_Filtered.fas', dicoFastapG4)

def build_arg_parser():
	GITDIR = os.getcwd()+'/'
	parser = argparse.ArgumentParser(description = 'FilterG4DependdingOnTr')
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	main(path)
