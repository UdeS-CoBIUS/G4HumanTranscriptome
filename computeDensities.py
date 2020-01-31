#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	This script will read results files from Wt and shuffled data set and
        compute all densities at a location by subclass level. The output of
        this script is a tsv file with all statistic computed : densities,
        percent of location with at least one pG4r, content of each nucleotide,
        GC content.

Command Line:
    For one chr : python ~/PATH/getMainDensities.py
"""

import os
import re
import argparse
import countLocation
import pandas as pd
from Bio import SeqIO
from pprint import pprint
import recurentFunction as rF

def importData(filename):
    """Creates from scratch a data frame with all nucleotides count for each Bt location.

    Creates an empty data frame and fill it by counting, for each location in
    the inputfile, the number of nucleotides.
    The input file contains all fasta sequences of all locations possible for
    one chromosome. For each locations, transcripts and biotypes are given in
    the id. This allows to get one row for each biotype locations.

    :param filename: fasta file of all possible locatinos for one chromosome.
	:type filename: string

    :returns: df, contains all biotype locations for one chromosome with their
        nulceotides content.Biotypes location are not yet unique.
    :rtype: dataFrame
    """
    df = pd.DataFrame(columns = ['LocID', 'Location', 'Biotype', 'nuclA', 'nuclT',
        'nuclG', 'nuclC', 'nuclN', 'nbTr'])
    dicoTmp = {}
    fastaOrigin = SeqIO.parse(open(filename),'fasta')
    for fasta in fastaOrigin:
        name, seq = fasta.id, str(fasta.seq)
        if name.split(':')[5]:
            location = name.split(':')[1]
            listTrBt = name.split(':')[5].split(';')[0].split('|')
            dicoTrBt = { TrBt.split('-')[0] : TrBt.split('-')[1] for TrBt in listTrBt}
            for tr in dicoTrBt:
                if not ((location == '3UTR' or location == '5UTR') and
                    rF.addTypeTr(dicoTrBt[tr]) != 'Coding'):
                    #if the annotation is good
                    LocID = location+'-'+dicoTrBt[tr]
                    if LocID not in dicoTmp:
                        dicoTmp[LocID] = {'LocID' : LocID,
                                        'Location' : location,
                                        'Biotype' : dicoTrBt[tr],
                                        'nuclA' : 0, 'nuclT' : 0,
                                        'nuclG' : 0, 'nuclC' : 0,
                                        'nuclN' : 0, 'nbTr' : [tr]}
                    dicoTmp[LocID].update({'nuclA' : dicoTmp[LocID]['nuclA'] + seq.count('A'),
                        'nuclT' : dicoTmp[LocID]['nuclT'] + seq.count('T'),
                        'nuclG' : dicoTmp[LocID]['nuclG'] + seq.count('G'),
                        'nuclC' : dicoTmp[LocID]['nuclC'] + seq.count('C'),
                        'nuclN' : dicoTmp[LocID]['nuclN'] + seq.count('N')})
                    dicoTmp[LocID]['nbTr'].append(tr)
    listTodf = []
    for locID in dicoTmp:
        listTodf.append(dicoTmp[locID])
    dfTmp = pd.DataFrame(listTodf)
    df = df.append(dfTmp)
    return(df)

def sumNt(df):
    """Sums all nucleotides content to get one row for each biotype location.

    The input dataFrame does not possess unique biotype location. To make them
    unique, all row of the same biotype location are summed together.

    :param df: contains all non unique biotype locations with there nucleotides
        content.
	:type df: dataFrame

    :returns: df, with unique row for eahc biotype locaions.
    :rtype: dataFrame
    """
    dfUniq = pd.DataFrame(columns = ['LocID', 'Location', 'Biotype', 'nuclA', 'nuclT',
        'nuclG', 'nuclC', 'nuclN', 'Tot', 'GC', 'nbTr'])
    groups = df.groupby('LocID')
    for name, group in groups:
        A = sum(list(group.nuclA))
        T = sum(list(group.nuclT))
        G = sum(list(group.nuclG))
        C = sum(list(group.nuclC))
        N = sum(list(group.nuclN))
        Tot = [ list(set(listTr))for listTr in group.nbTr ]
        Tot = [ len(liste) for liste in Tot]
        Tot = sum(Tot)
        row = {'LocID' : name,
                'Location' : list(set(group.Location)),
                'Biotype' : list(set(group.Biotype)),
                'nuclA' : A, 'nuclT' : T, 'nuclG' : G, 'nuclC' : C,
                'nuclN' : N, 'Tot' : A+T+G+C+N,
                'GC' : float(G+C)/(A+T+G+C+N)*100,
                'nbTr' : Tot}
        row = pd.DataFrame(row, index=[len(dfUniq)+1])
        dfUniq = dfUniq.append(row)
    return(dfUniq)

def getGCFromFile(path):
    """Creates from scratch a data frame with all nucleotides count for each Bt.

    Creates an empty data frame and fill it after reader each fasta shuffled
    file. The data frame contains one columns for each nucleotides.

    :param path: path to all file needed.
	:type path: string

    :returns: df, containing all nucleotides count for each subclasses locations.
    :rtype: dataFrame
    """
    df = pd.DataFrame(columns = ['LocID', 'Location', 'Biotype', 'nuclA',
        'nuclT', 'nuclG', 'nuclC', 'nuclN'])
    for path2, dirs, files in os.walk(path):
        for chr in dirs:
            if 'chr' in chr:
                filename = path2+chr+'/Shuffle_'+chr+'_location.fas'
                df = df.append(importData(filename))
                filename = path2+chr+'/Shuffle_'+chr+'_junction.fas'
                df = df.append(importData(filename))
    df = sumNt(df)
    return(df)

def getDicoNbpG4rShuf(pG4rFile):
    """Import from a file the pG4r number and the number of transcript with it.

    This function is made for the shuffled dataset.
    This function aims to import into a dictionary the number of pG4r for each
    biotype locations, but also the number of transcript with at least one pG4r.
    To do that, we check each row of the file and then add to the dictionary
    this pG4r ( +1 for the number of pG4r and +Id for the number of transcript).
    For the nmber of transcript the list with all ID is make unique and then its
    length gave the number of transcript with at least one pG4r.

    :param pG4rFile: pG4r filename for the shuffled dataset.
	:type pG4rFile: string

    :returns: dicopG4r, {NbG4 : {LocID : nbpG4r}, nbTrWithpG4 : {LocID : nbTr}}.
    :rtype: dictionary
    """
    dicopG4r = {'NbG4' : {},
                'nbTrWithpG4' : {}}
    with open(pG4rFile) as f:
        lines = f.read().splitlines()
        for l in lines:
            l = l.rstrip()
            words = l.split('\t')
            if words[0] != 'pG4rID' and words[0]:
                id = words[0].split(';')[0]
                if id.split(':')[5]:
                    location = id.split(':')[1]
                    listTrBt = id.split(':')[5].split(';')[0].split('|')
                    listBt = [ TrBt.split('-')[1] for TrBt in listTrBt ]
                    dicoBt = { bt : [] for bt in listBt}
                    for TrBt in listTrBt:
                        dicoBt[ TrBt.split('-')[1] ].append(TrBt.split('-')[0])
                    for Bt in listBt:
                        locID = location+'-'+Bt
                        if locID not in dicopG4r['NbG4']:
                            dicopG4r['NbG4'][locID] = 0
                        if locID not in dicopG4r['nbTrWithpG4']:
                            dicopG4r['nbTrWithpG4'][locID] = []
                        dicopG4r['NbG4'][locID] += 1
                        dicopG4r['nbTrWithpG4'][locID].extend(dicoBt[Bt])
    for locID in dicopG4r['nbTrWithpG4']:
        dicopG4r['nbTrWithpG4'][locID] = len(list(set(dicopG4r['nbTrWithpG4'][locID])))
    return(dicopG4r)

def changeLocName(location):
    """Changes names of locations.

    This function aims to unified names location between the shuffled dataset
    and the wildtype dataset.

    :param location: wt name location.
	:type location: string

    :returns: location, updated with the 'good' name.
    :rtype: string
    """
    if location == 'junction_IntronNC_ExonNC' or location == 'junction_Intron_CDS':
        location = 'acceptor'
    elif location == 'junction_ExonNC_IntronNC' or location == 'junction_CDS_Intron':
        location = 'donor'
    elif location == 'junction_CDS_CDS' or location == 'junction_ExonNC_ExonNC':
        location = 'junction'
    elif location == '3':
        location = '3UTR'
    elif location == '5':
        location = '5UTR'
    elif location == 'ExonNC':
        location = 'exon'
    elif location == 'Intron' or location == 'IntronNC':
        location = 'intron'
    elif location == 'junction_5_CDS':
        location = 'StartCodon'
    elif location == 'junction_CDS_3':
        location = 'StopCodon'
    return location

def getDicoNbpG4rWt(pG4rFile):
    """Import from a file the pG4r number and the number of transcript with it.

    This function is made for the Wt dataset. Exon location are added because
    they are not explicit (5UTR + CDS + 3UTR + codons). Also, location names
    are changed to get the one used in the analysis.
    This function aims to import into a dictionary the number of pG4r for each
    biotype locations, but also the number of transcript with at least one pG4r.
    To do that, we check each row of the file and then add to the dictionary
    this pG4r ( +1 for the number of pG4r and +Id for the number of transcript).
    For the nmber of transcript the list with all ID is make unique and then its
    length gave the number of transcript with at least one pG4r.

    :param pG4rFile: pG4r filename for the Wt dataset.
	:type pG4rFile: string

    :returns: dicopG4r, {NbG4 : {LocID : nbpG4r}, nbTrWithpG4 : {LocID : nbTr}}.
    :rtype: dictionary
    """
    dicopG4r = {'NbG4' : {},
                'nbTrWithpG4' : {}}
    exonBtMissing = ['nonsense_mediated_decay', 'protein_coding',
                    'IG_V_gene', 'non_stop_decay']
    with open(pG4rFile) as f:
        lines = f.read().splitlines()
        for l in lines:
            l = l.rstrip()
            words = l.split('\t')
            if words[0] != 'pG4rID' and words[0]:
                id = words[0]
                location = words[5]
                Bt = words[6]
                tr = id.split('|')[0]
                location = changeLocName(location)
                locID = location+'-'+Bt
                if locID not in dicopG4r['NbG4']:
                    dicopG4r['NbG4'][locID] = 0
                if locID not in dicopG4r['nbTrWithpG4']:
                    dicopG4r['nbTrWithpG4'][locID] = []
                dicopG4r['NbG4'][locID] += 1
                dicopG4r['nbTrWithpG4'][locID].append(tr)
                if location in ['5UTR', 'CDS', '3UTR', 'StartCodon', 'StopCodon'] and Bt in exonBtMissing:
                    locID = 'exon-'+Bt
                    if locID not in dicopG4r['NbG4']:
                        dicopG4r['NbG4'][locID] = 0
                    if locID not in dicopG4r['nbTrWithpG4']:
                        dicopG4r['nbTrWithpG4'][locID] = []
                    dicopG4r['NbG4'][locID] += 1
                    dicopG4r['nbTrWithpG4'][locID].append(tr)
    for locID in dicopG4r['nbTrWithpG4']:
        dicopG4r['nbTrWithpG4'][locID] = len(list(set(dicopG4r['nbTrWithpG4'][locID])))
    return(dicopG4r)

def addNbpG4r(dico, df, type):
    """Adds pG4r number for one dataset.

    The input dataFrame is browsed, and pG4r number and transcript number are
    added one by one for each biotype location.

    :param dico: contains the pG4r number for one dataset.
	:type dico: dictionary
    :param df: contains all biotype locations and their nucleotides content.
	:type df: dataFrame
    :param type: type of the dataset (Wt or shuffled).
	:type type: string

    :returns: df, updated with the number of pG4r for one dataset.
    :rtype: dataFrame
    """
    df = df.reset_index()
    name = 'NbpG4r'+type
    name2 = 'NbTrpG4'+type
    df[name] = 0
    df[name2] = 0
    for index, row in df.iterrows():
        if row.LocID in dico['NbG4']:
            df[name].iloc[index] = dico['NbG4'][row.LocID]
        if row.LocID in dico['nbTrWithpG4']:
            df[name2].iloc[index] = dico['nbTrWithpG4'][row.LocID]
    return df

def addpG4rNumber(pG4rFileShuf, pG4rFileWt, df):
    """Adds to the input dataFrame the number of pG4r for the wt and shuffles DS.

    First import information about pG4r (parsed file of the G4RNA screener
    ouput file). Then, those information are added by calling addNbpG4r for
    each dataset.

    :param pG4rFileShuf: pG4r Shuffled filename.
	:type pG4rFileShuf: string
    :param pG4rFileWt: pG4r Wt filename.
	:type pG4rFileWt: string
    :param df: contains all biotype locations and their nucleotides content.
	:type df: dataFrame

    :returns: df, updated with the number of pG4r for each data set.
    :rtype: dataFrame
    """
    dicoShuf = getDicoNbpG4rShuf(pG4rFileShuf)
    dicoWt = getDicoNbpG4rWt(pG4rFileWt)
    df = addNbpG4r(dicoShuf, df, 'Shuf')
    df = addNbpG4r(dicoWt, df, 'Wt')
    return df

def computeDensities(df, type):
    """Computes pG4r densities.

    Browses the input dataFrame, and for each row computes the pG4r density.
    If the row is a segment location it is : numberOfpG4r / lengthTot * 1000,
    and if it is a point location : numberOfpG4r / numberTotOfLocation.

    :param df: contains nucleotides content, number of pG4r and number of
        number transcripts and number of locations.
	:type df: dataFrame
    :param type: name of the dataset (Wt or Shuf)
	:type type:string

    :returns:df, updated with densities.
    :rtype: dataFrame
    """
    namepG4r = 'NbpG4r'+type
    nameDensity = 'Density'+type
    df[nameDensity] = 0.0
    for index, row in df.iterrows():
        if row[namepG4r] != 0:
            if row.Type == 'Segment':
                d = float(row[namepG4r])/row.Tot * 1000
            else:
                d = float(row[namepG4r])/row.NbLocation
        else:
            d = 0
        df[nameDensity].iloc[index] = d
    return df

def addTypeTr(biotype):
    """Apply function to retrieve the type of a transcript depending on biotype.

    :param biotype: biotype of a transcript.
    :type biotype: string

    :returns: type, coding or non coding.
    :rtype: string
    """
    dicoFam = rF.createDicoFamily()
    for bt in dicoFam:
        if biotype in dicoFam[bt]:
            return bt

def addType(location):
    """Returns the location type, made for an apply.

    :param location:
	:type location: string

    :returns: type of the location.
    :rtype: string
    """
    if location in ['donor', 'acceptor', 'junction', 'StartCodon', 'StopCodon']:
        return 'Point'
    else:
        return 'Segment'

def addNbLocation(path, df):
    """Adds the number of locations with at least one pG4r in it.

    First all number are retrieved : number tot of location by biotype, number
    of location with at least one pG4r for each dataset. Then they are mapped
    to the dataFrame with nucleotides content and the number of pG4r.

    :param path: path to all needed file.
	:type path: string
    :param df: contains nucleotides content, number of pG4r and number of
        number transcripts.
	:type df: dataFrame

    :returns: df, updated with the numbers of location with or without pG4r for
        each dataset.
    :rtype: dataFrame
    """
    tmp = pd.DataFrame()
    dicoTot = countLocation.getAllLocNumMulti(path, 'Bt')
    dicoWt = countLocation.getpG4WtLocNum(path, 'Bt')
    dicoShuf = countLocation.getpG4ShufLocNum(path, 'Bt')
    df['NbLocation'] = df['LocID'].map(dicoTot)
    df['NbpG4rLocWt'] = df['LocID'].map(dicoWt)
    df['NbpG4rLocShuf'] = df['LocID'].map(dicoShuf)
    return df

def main(path):
    """ Main function to obtain all location-subclasses densities.

    First a dataframe with all nucleotides count for each sublclasses is created.
    The the pG4r is added for Wt and shuffled data set. Then, we add the
    number of unique location that contain at least one pG4r. Finally, densities
    are computed.
    This function return nothing but create a file : TotDataDensities.csv.

    :param path: path to all file needed.
	:type path: string
    """
    pG4rFileShuf = path+'Resulst/All/pG4r_shuffle.csv'
    pG4rFileWt = path+'Results/All/HS_All_G4InTranscript.txt'
    df = getGCFromFile(path+'Data/')
    df = addpG4rNumber(pG4rFileShuf, pG4rFileWt, df)
    df['Class'] = df.Biotype.apply(addTypeTr)
    df['Type'] = df.Location.apply(addType)
    df = addNbLocation(path, df)
    df = computeDensities(df, 'Shuf')
    df = computeDensities(df, 'Wt')
    df = df.fillna(0)
    del df['level_0']
    df = df.drop_duplicates(subset=None, keep='first', inplace=False)
    df.to_csv(path_or_buf=path+'/Results/All/TotDataDensitesTest.csv', header=True, index=None, sep='\t')

def build_arg_parser():
    parser = argparse.ArgumentParser(description = 'analyseGC')
    GITDIR = os.getcwd()+'/'
    parser.add_argument ('-p', '--path', default = GITDIR)
    return parser

if __name__ == '__main__':
    parser = build_arg_parser()
    arg = parser.parse_args()
    path = arg.path
    main(path)
