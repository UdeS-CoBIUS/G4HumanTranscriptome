#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	August 2019

Description:
	This file is more like a library called when statistics are computed (
    getMainDensities and getDataFig). Its aims is to compute the number of
    transcript with at least one pG4r in it. A special script is made for
    because we need to compute each time, otherwise we could get the same
    transcript in different condition and count it 2 times.
"""

import os
import argparse
from pprint import pprint
import recurentFunction as rF

def importTr(filename):
    """Imports in a list, all transcripts id. IDs are uniq.

    :param filename: name of the file containing all transcripts.
	:type filename: string

    :returns: listTr, a list of uniq transcripts id.
    :rtype: list
    """
    listTr = []
    with open(filename) as f:
        lines = f.read().splitlines()
        for l in lines:
            l = l.rstrip()
            words = l.split('|')
            listTr.append(words[0])
    return (list(set(listTr)))

def readShufpG4r(pG4rFile):
    """Reads a file containing all pG4r Shuf and get their transcripts.

    This function aims to retrieve transcript id of those that have a G4
    predicted in them.

    :param filename: name of the file containing all shuf pG4r.
	:type filename: string

    :returns: listTr, contains all uniq tr id that have at least one pG4r.
    :rtype: list
    """
    listTr = []
    with open(pG4rFile) as f:
        lines = f.read().splitlines()
        for l in lines:
            l = l.rstrip()
            words = l.split('\t')
            if words[0] != 'pG4rID' and words[0]:
                id = words[0].split(';')[0]
                if id.split(':')[5]:
                    if id.split(':')[1] in ['exon', 'intron']:
                        listTrBt = id.split(':')[5].split(';')[0].split('|')
                        listTr.extend( [ TrBt.split('-')[0] for TrBt in listTrBt ] )
    return list(set(listTr))

def readWtpG4r(filename):
    """Reads a file containing all pG4r WT and get their transcripts.

    This function aims to retrieve transcript id of those that have a G4
    predicted in them.

    :param filename: name of the file containing all wt pG4r.
	:type filename: string

    :returns: listTr, contains all uniq tr id that have at least one pG4r.
    :rtype: list
    """
    listTr = []
    with open(filename) as f:
        lines = f.read().splitlines()
        for l in lines:
            l = l.rstrip()
            words = l.split('|')
            listTr.append(words[0])
    return list(set(listTr))

def getNumberTrinClass(dicoBt, listTr):
    """Gets the number of transcripts in all classes.

    This function aims to first print the number of transcripts in all classes,
    and second to return a dictionary containing all transcript and their class.

    :param dicoBt: contains all transcript and their biotype {idTr : biotype}.
	:type dicoBt: dictionary
    :param listTr: contains all uniq transcripts id, that have a biotype in the
        human transcriptome.
	:type listTr: list

    :returns: dicoTrClass, dictionary with transcipt and their class {idTr : class}.
    :rtype: dictionary
    """
    dicoClass = rF.createDicoFamilyFiltered()
    dicoRes = {'Coding':0,'Pseudogene':0,'LongNC':0,'ShortNC':0,'Predictif':0}
    dicoTrClass = {}
    for tr in listTr:
        for classe in dicoClass:
            if dicoBt[tr] in dicoClass[classe]:
                dicoRes[classe] += 1
                if tr not in dicoTrClass:
                    dicoTrClass[tr] = classe
    return dicoTrClass, dicoRes

def getNumberTrinBt(dicoBt, listTr):
    """Gets the number of transcripts in all subclasses.

    :param dicoBt: contains all transcript and their biotype {idTr : biotype}.
	:type dicoBt: dictionary
    :param listTr: contains all uniq transcripts id, that have a biotype in the
        human transcriptome.
	:type listTr: list
    """
    dicoClass = rF.createDicoFamilyFiltered()
    dicoRes = {}
    dicoTrBt = {}
    for tr in listTr:
        for classe in dicoClass:
            if dicoBt[tr] in dicoClass[classe]:
                if dicoBt[tr] not in dicoRes:
                    dicoRes[ dicoBt[tr] ] = 0
                dicoRes[ dicoBt[tr] ] += 1
    return dicoRes

def countpG4rByClass(dicoTrClass, listTrpG4r):
    """Counts the number of transcripts with pG4r for each transcript class.

    :param dicoTrClass: dictionary with transcipt and their class {idTr : class}.
	:type dicoTrClass: dictionary
    :param listTrpG4r: contains all uniq tr id that have at least one pG4r.
	:type listTrpG4r: list
    """
    dicoRes = {'Coding':0,'Pseudogene':0,'LongNC':0,'ShortNC':0,'Predictif':0}
    for tr in listTrpG4r:
        if tr in dicoTrClass:
                dicoRes[ dicoTrClass[tr] ] += 1
    return dicoRes

def countpG4rByBt(dicoTrBt, listTrpG4r):
    """Counts the number of transcripts with pG4r for each transcript subclass.

    :param dicoTrBt: dictionary with transcipt and their subclass {idTr : subclass}.
	:type dicoTrBt: dictionary
    :param listTrpG4r: contains all uniq tr id that have at least one pG4r.
	:type listTrpG4r: list
    """
    dicoRes = {}
    for tr in listTrpG4r:
        if tr in dicoTrBt:
            if dicoTrBt[tr] not in dicoRes:
                dicoRes[ dicoTrBt[tr] ] = 0
            dicoRes[ dicoTrBt[tr] ] += 1
    return dicoRes

def getFig5Percent(path):
    fileTr = path+'Data/transcriptType/HS_transcript_unspliced_All.txt'
    fileBt = path+'Data/transcriptType/transcriptType_All.txt'
    filepG4rWt = path+'Results/All/HS_All_G4InTranscript.txt'
    filepG4rShuf = path+'Results/All/pG4r_shuffle.csv'
    dicoNbTr = { 'Wt' : {}, 'Shuf' : {}, 'Tot' : {} }

    dicoBt = rF.createDictionaryBiotypeByTranscript(fileBt)
    listTr = importTr(fileTr)
    # we keep only transcript with a biotype known.
    listTot = [ tr for tr in listTr if tr in dicoBt]
    dicoNbTrClass = getNumberTrinBt(dicoBt, listTot)
    dicoNbTr['Tot'].update(dicoNbTrClass)

    listTrpG4rWt = readWtpG4r(filepG4rWt)
    dicoNbTrWt = countpG4rByBt(dicoBt, listTrpG4rWt)
    dicoNbTr['Wt'].update(dicoNbTrWt)

    listTrpG4rShuf = readShufpG4r(filepG4rShuf)
    dicoNbTrShuf = countpG4rByBt(dicoBt, listTrpG4rShuf)
    dicoNbTr['Shuf'].update(dicoNbTrShuf)
    return dicoNbTr

def getFig3Percent(path):
    fileTr = path+'Data/transcriptType/HS_transcript_unspliced_All.txt'
    fileBt = path+'Data/transcriptType/transcriptType_All.txt'
    filepG4rWt = path+'Results/All/HS_All_G4InTranscript.txt'
    filepG4rShuf = path+'Results/All/pG4r_shuffle.csv'
    dicoNbTr = { 'Wt' : {}, 'Shuf' : {}, 'Tot' : {} }

    dicoBt = rF.createDictionaryBiotypeByTranscript(fileBt)
    listTr = importTr(fileTr)
    # we keep only transcript with a biotype known.
    listTot = [ tr for tr in listTr if tr in dicoBt]
    dicoTrClass, dicoNbTrClass = getNumberTrinClass(dicoBt, listTot)
    dicoNbTr['Tot']['Global'] = len(set(listTot))
    dicoNbTr['Tot'].update(dicoNbTrClass)

    listTrpG4rWt = readWtpG4r(filepG4rWt)
    dicoNbTrWt = countpG4rByClass(dicoTrClass, listTrpG4rWt)
    dicoNbTr['Wt']['Global'] = len(listTrpG4rWt)
    dicoNbTr['Wt'].update(dicoNbTrWt)

    listTrpG4rShuf = readShufpG4r(filepG4rShuf)
    dicoNbTrShuf = countpG4rByClass(dicoTrClass, listTrpG4rShuf)
    dicoNbTr['Shuf']['Global'] = len(listTrpG4rShuf)
    dicoNbTr['Shuf'].update(dicoNbTrShuf)
    return dicoNbTr

def getSupPercent(path):
    fileTr = path+'Data/transcriptType/HS_transcript_unspliced_All.txt'
    fileBt = path+'Data/transcriptType/transcriptType_All.txt'
    filepG4rWt = path+'Results/All/HS_All_G4InTranscript.txt'
    filepG4rShuf = path+'Results/All/pG4r_shuffle.csv'
    dicoNbTr = { 'Wt' : {}, 'Shuf' : {}, 'Tot' : {} }

    dicoBt = rF.createDictionaryBiotypeByTranscript(fileBt)
    listTr = importTr(fileTr)
    # we keep only transcript with a biotype known.
    listTot = [ tr for tr in listTr if tr in dicoBt]
    dicoTrBt = getNumberTrinBt(dicoBt, listTot)
    dicoNbTr['Tot'].update(dicoTrBt)

    listTrpG4rWt = readWtpG4r(filepG4rWt)
    dicoNbTrWt = countpG4rByBt(dicoBt, listTrpG4rWt)
    dicoNbTr['Wt'].update(dicoNbTrWt)

    listTrpG4rShuf = readShufpG4r(filepG4rShuf)
    dicoNbTrShuf = countpG4rByBt(dicoBt, listTrpG4rShuf)
    dicoNbTr['Shuf'].update(dicoNbTrShuf)
    return dicoNbTr

def main(path):
    # fileTr = path+'transcriptType/HS_transcript_unspliced_All.txt'
    # fileBt = path+'transcriptType/transcriptType_All.txt'
    # filepG4rWt = path+'chrAll/HS_All_G4InTranscript.txt'
    # filepG4rShuf = path+'chrAll/pG4r_shuffle5.csv'
    # dicoBt = rF.createDictionaryBiotypeByTranscript(fileBt)
    # listTr = importTr(fileTr)
    # # we keep only transcript with a biotype known.
    # listTot = [ tr for tr in listTr if tr in dicoBt]
    # print('Number tot of transcripts : ', len(set(listTot)))
    #
    # dicotrClass = getNumberTrinClass(dicoBt, listTot)
    # listTrpG4rWt = readWtpG4r(filepG4rWt)
    # print('Number of transcripts with WT pG4r: ', len(listTrpG4rWt))
    # countpG4rByClass(dicotrClass, listTrpG4rWt, 'WT')
    # listTrpG4rShuf = readShufpG4r(filepG4rShuf)
    # print('Number of transcripts with Shuf pG4r: ', len(listTrpG4rShuf))
    # countpG4rByClass(dicotrClass, listTrpG4rShuf, 'Shuf')
    #
    # getNumberTrinBt(dicoBt, listTot)
    # countpG4rByBt(dicoBt, listTrpG4rWt, 'WT')
    # countpG4rByBt(dicoBt, listTrpG4rShuf, 'Shuf')
    # getSupPercent(path)
    pprint(getFig3Percent(path))

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'analyseGC')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    mainPath = arg.path
    main(mainPath)
