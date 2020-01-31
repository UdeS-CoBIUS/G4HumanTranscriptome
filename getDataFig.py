#!/usr/bin/env python
# -*- coding: utf-8 -*-:

"""

Copyright:
	Copyright Universite of Sherbrooke, departement of biochemistry and
	departement	of computation.

Date:
	September 2019

Description:
	With the output file of getMainDensities.py, this script compute all other
		statistics needed. There is one tsv output by figure.

Command Line:
	For one chr : python ~/PATH/getDataFig.py
"""

import os
import argparse
import pandas as pd
import countTranscript
import countLocation
from pprint import pprint
import recurentFunction as rF

def removeBiotype(df):
	"""Removes from a data frame, informations of some biotype we don't want.

	A lot of biotypes are removed for mainly 2 reasons : the biotype doesn't
	have enough transcript (less than 50) or the biotype by definition is linked
	to a high DNA/RNA variation (like IG genes).

	:param df: origin dataframe with all biotypes.
	:type df: DataFrame

	:returns: updated dataframe with only biotype of interest
	:rtype: DataFrame
	"""
	df = df[ df.Biotype != 'IG_C_gene']
	df = df[ df.Biotype != 'IG_D_gene']
	df = df[ df.Biotype != 'IG_J_gene']
	df = df[ df.Biotype != 'IG_V_gene']
	df = df[ df.Biotype != 'pseudogene']
	df = df[ df.Biotype != 'rRNA']
	df = df[ df.Biotype != 'sRNA']
	df = df[ df.Biotype != 'TR_C_gene']
	df = df[ df.Biotype != 'TR_D_gene']
	df = df[ df.Biotype != 'TR_J_gene']
	df = df[ df.Biotype != 'TR_V_gene']
	df = df[ df.Biotype != 'macro_lncRNA']
	df = df[ df.Biotype != 'bidirectional_promoter_lncRNA']
	df = df[ df.Biotype != '3prime_overlapping_ncRNA']
	df = df[ df.Biotype != 'non_coding']
	df = df[ df.Biotype != 'pseudogene']
	df = df[ df.Biotype != 'TR_J_pseudogene']
	df = df[ df.Biotype != 'IG_C_pseudogene']
	df = df[ df.Biotype != 'IG_J_pseudogene']
	df = df[ df.Biotype != 'IG_pseudogene']
	df = df[ df.Biotype != 'TR_V_pseudogene']
	df = df[ df.Biotype != 'polymorphic_pseudogene']
	df = df[ df.Biotype != 'IG_V_pseudogene']
	df = df[ df.Biotype != 'TEC']
	df = df[ df.Biotype != 'Predictif']
	df = df[ df.Biotype != 'ribozyme']
	df = df[ df.Biotype != 'scRNA']
	df = df[ df.Biotype != 'scaRNA']
	df = df[ df.Biotype != 'snRNA']
	df = df[ df.Biotype != 'snoRNA']
	df = df[ df.Biotype != 'vaultRNA']
	df = df[ df.Biotype != 'translated_processed_pseudogene']
	return df

def sumSubTable(group, name):
	NbpG4rShuf = sum(group.NbpG4rShuf)
	NbpG4rWt = sum(group.NbpG4rWt)
	NbTrpG4Wt = sum(group.NbTrpG4Wt)
	NbTrpG4Shuf = sum(group.NbTrpG4Shuf)
	nbTr = sum(group.nbTr)
	Tot = sum(group.Tot)
	G = sum(group.nuclG)
	C = sum(group.nuclC)
	DensityShuf = float(NbpG4rShuf)/Tot *1000
	DensityWt = float(NbpG4rWt)/Tot *1000
	row = {'nuclG' : G,
			'nuclC' : C,
			'nbTr' : nbTr,
			'NbpG4rWt' : NbpG4rWt,
			'NbpG4rShuf' : NbpG4rShuf,
			'NbTrpG4Wt' : NbTrpG4Wt,
			'NbTrpG4Shuf' : NbTrpG4Shuf,
			'Tot' : Tot}
	return row

def computePercentage(df, type):
	tmp = pd.DataFrame()
	tmp = tmp.append(df)
	locType = 'NbpG4rLoc'+type
	namePercent = 'Percent'+type
	tmp[namePercent] = 0.0
	for index, row in tmp.iterrows():
		locId = row.LocID
		if row[locType] != 0:
				d = float(row[locType])/row.NbLocation *100
		else:
			d = 0
		tmp.loc[ tmp['LocID'] == locId, [namePercent] ] = d
		# tmp[namePercent].iloc[index] = d
	return tmp

def computeDensity(df, type):
	tmp = pd.DataFrame()
	tmp = tmp.append(df)
	tmp['GC'] = (tmp['nuclG'] + tmp['nuclC']) / tmp['Tot'] *100
	if type == 'Point':
		tmp['DensityWt'] = tmp['NbpG4rWt'] / tmp['NbLocation']
		tmp['DensityShuf'] = tmp['NbpG4rShuf'] / tmp['NbLocation']
	else:
		tmp['DensityWt'] = tmp['NbpG4rWt'] / tmp['Tot'] * 1000
		tmp['DensityShuf'] = tmp['NbpG4rShuf'] / tmp['Tot'] * 1000
	return tmp

def computePercent(dico):
	if dico['nbTr'] == 0:
		dico['PercentWt'] = 0
		dico['PercentShuf'] = 0
	else:
		dico['PercentWt'] = float(dico['NbTrpG4Wt']) / float(dico['nbTr']) * 100
		dico['PercentShuf'] = float(dico['NbTrpG4Shuf']) / float(dico['nbTr']) * 100
	return dico

def getFigBysubClass(df, path, nameClass):
	"""Same as for figure 3 but for subclass level.

	This function aims to filter the input dataFrame to get information about
	pG4r at class and subclass level. The final dataFrame contains
	information to make lolipop plot for GC content, densities and number
	of transcript.
	The glaobal statistic are cuomputed for each location and each informations
	(GC, densities, percent).

	:param df: contains all informations for all biotype locations.
	:type df: dataFrame
	:param path: general path with all files.
	:type path: string

	:returns: classDf, contains data for a location figure.
	:rtype: dataFrame
	"""
	tmp = pd.DataFrame()
	tmp = tmp.append(df)
	dicoNbTrClass = countTranscript.getFig3Percent(path)
	dicoNbTrBt = countTranscript.getFig5Percent(path)
	del tmp['nuclA']
	del tmp['nuclT']
	del tmp['nuclN']
	del tmp['Type']
	classDf = pd.DataFrame()
	classDftmp = tmp[ tmp.Class == nameClass]
	groups = classDftmp.groupby('Biotype')
	for name, group in groups:
		groupFilter = group[ group.Location == 'intron' ]
		groupFilter = groupFilter.append( group[ group.Location == 'exon' ])
		row = sumSubTable(groupFilter, name)
		row['Biotype'] = name
		row['Class'] = nameClass
		if name not in dicoNbTrBt['Tot']:
			dicoNbTrBt['Tot'][name] = 0
		if name not in dicoNbTrBt['Wt']:
			dicoNbTrBt['Wt'][name] = 0
		if name not in dicoNbTrBt['Shuf']:
			dicoNbTrBt['Shuf'][name] = 0
		row['nbTr'] = dicoNbTrBt['Tot'][name]
		row['NbTrpG4Wt'] = dicoNbTrBt['Wt'][name]
		row['NbTrpG4Shuf'] = dicoNbTrBt['Shuf'][name]
		row.update(computePercent(row))
		row = pd.DataFrame(row, index=[len(classDftmp)+1])
		classDf = classDf.append(row)
	row = {'Class' : nameClass,
			'Biotype' : nameClass,
			'nuclG' : sum(classDftmp.nuclG),
			'nuclC' : sum(classDftmp.nuclC),
			'nbTr' : dicoNbTrClass['Tot'][nameClass],
			'NbpG4rWt' : sum(classDftmp.NbpG4rWt),
			'NbpG4rShuf' : sum(classDftmp.NbpG4rShuf),
			'NbTrpG4Wt' : dicoNbTrClass['Wt'][nameClass],
			'NbTrpG4Shuf' : dicoNbTrClass['Shuf'][nameClass],
			'Tot' : sum(classDftmp.Tot)}
	row.update(computePercent(row))
	row = pd.DataFrame(row, index=[len(classDf)+1])
	classDf = classDf.append(row)
	classDf = computeDensity(classDf, 'Segment')
	return classDf

def getFigDensityByLocation(df, path, nameClass):
	"""Gets the data frame for figure with pG4r densities by location.

	This function aims to filter the input dataFrame to get information about
	pG4r at location level. The final dataFrame contains information to make
	lolipop plot for GC content, densities and number of location figure for
	segmental and point locations
	The glaobal statistic are cuomputed for each location and each informations
	(GC, densities, percent).

	:param df: contains all informations for all biotype locations.
	:type df: dataFrame
	:param path: general path with all files.
	:type path: string
	:param nameClass: name of the class we want to get the data (Coding, LongNC,
		Pseudogene).
	:type nameClass: string

	:returns: classDf, contains data for a location figure.
	:rtype: dataFrame
	"""
	tmp =  pd.DataFrame()
	dicoTotMulti = countLocation.getAllLocNumMulti(path, 'Class')
	dicoWt = countLocation.getpG4WtLocNum(path, 'Class')
	dicoShuf = countLocation.getpG4ShufLocNum(path, 'Class')
	classDf = df[ df.Class == nameClass]
	if nameClass in ['Pseudogene', 'LongNC'] :
		#Those class should not get those location
		classDf = classDf[ classDf.Location != 'StartCodon']
		classDf = classDf[ classDf.Location != 'StopCodon']
	del classDf['nuclA']
	del classDf['nuclT']
	del classDf['nuclN']
	del classDf['Class']
	groups = classDf.groupby('Location')
	for name, group in groups:
		row = sumSubTable(group, name)
		row['Biotype'] = nameClass
		row['LocID'] = name+'-'+nameClass
		row['Location'] = name
		row['NbLocation'] = dicoTotMulti[name+'-'+nameClass]
		row['NbpG4rLocWt'] = dicoWt[name+'-'+nameClass]
		row['NbpG4rLocShuf'] = dicoShuf[name+'-'+nameClass]
		row['PercentWt'] = float(dicoWt[name+'-'+nameClass]) / float(dicoTotMulti[name+'-'+nameClass]) * 100
		row['PercentShuf'] = float(dicoShuf[name+'-'+nameClass]) / float(dicoTotMulti[name+'-'+nameClass]) * 100
		if name in ['donor', 'acceptor', 'junction', 'StartCodon', 'StopCodon']:
			row['Type'] = 'Point'
		else:
			row['Type'] = 'Segment'
		row = pd.DataFrame(row, index=[len(tmp)+1])
		tmp = tmp.append(row)
	tmp1 = computeDensity(tmp[ tmp.Type == 'Point'], 'Point')
	tmp2 = computeDensity(tmp[ tmp.Type == 'Segment'], 'Segment')
	classDf = classDf.append(tmp1)
	classDf = classDf.append(tmp2)
	classDf = computePercentage(classDf, 'Wt')
	classDf = computePercentage(classDf, 'Shuf')
	return classDf

def getFig3Data(df, path):
	"""Gets the data frame for figure 3.

	This function aims to filter the input dataFrame to get information about
	pG4r at global level and class level. The final dataFrame contains
	information to make lolipop plot for GC content, densities and number
	of transcript.
	The glaobal statistic are cuomputed for each location and each informations
	(GC, densities, percent).

	:param df: contains all informations for all biotype locations.
	:type df: dataFrame
	:param path: general path with all files.
	:type path: string

	:returns: classDf, contains data for a location figure.
	:rtype: dataFrame
	"""
	tmp = pd.DataFrame()
	# tmp = tmp.append(df)
	tmp = tmp.append(df[df.Location == 'exon'])
	tmp = tmp.append(df[df.Location == 'intron'])
	# print(df[df.Location == 'exon'].NbpG4rWt)
	# print(df[df.Location == 'intron'].NbpG4rWt)
	dicoNbTr = countTranscript.getFig3Percent(path)
	Global = pd.DataFrame()
	groups = tmp.groupby('Class')
	for name, group in groups:
		row = sumSubTable(group, name)
		row['Class'] = name
		row = pd.DataFrame(row, index=[len(Global)+1])
		Global = Global.append(row)
	# print(sum(Global.NbpG4rWt))
	row = {'Class' : 'Global',
			'nuclG' : sum(Global.nuclG),
			'nuclC' : sum(Global.nuclC),
			'NbpG4rWt' : sum(Global.NbpG4rWt),
			'NbpG4rShuf' : sum(Global.NbpG4rShuf),
			'Tot' : sum(Global.Tot)}
	row = pd.DataFrame(row, index=[len(Global)+1])
	Global = Global.append(row)
	Global['nbTr'] = Global['Class'].map( dicoNbTr['Tot'] )
	Global['NbTrpG4Wt'] = Global['Class'].map( dicoNbTr['Wt'] )
	Global['NbTrpG4Shuf'] = Global['Class'].map( dicoNbTr['Shuf'] )
	Global['PercentWt'] = Global['NbTrpG4Wt'] / Global['nbTr'] * 100
	Global['PercentShuf'] = Global['NbTrpG4Shuf'] / Global['nbTr'] * 100
	Global = computeDensity(Global, 'Segment')
	return Global

def main(path):
	results = path+'Results/All/TotDataDensites.csv'
	try:
		df = pd.read_csv(results, sep='\t', index_col=0)
	except:
		print("This file couldn't be converted in data frame : " + results)
	else:
		df = removeBiotype(df)
		df = df.reset_index()
		fig3 = getFig3Data(df, path)
		fig3.to_csv(path_or_buf=path+'/Results/Figures/DataFig3.csv', header=True, index=None, sep='\t')
		fig4 = getFigBysubClass(df, path, 'Coding')
		fig4.to_csv(path_or_buf=path+'/Results/Figures/DataFig4.csv', header=True, index=None, sep='\t')
		fig5 = getFigDensityByLocation(df, path, 'Coding')
		fig5.to_csv(path_or_buf=path+'/Results/Figures/DataFig5.csv', header=True, index=None, sep='\t')
		fig6 = getFigBysubClass(df, path, 'LongNC')
		fig6 = fig6.append(getFigBysubClass(df, path, 'ShortNC'))
		fig6.to_csv(path_or_buf=path+'/Results/Figures/DataFig6.csv', header=True, index=None, sep='\t')
		fig7 = getFigDensityByLocation(df, path, 'LongNC')
		fig7.to_csv(path_or_buf=path+'/Results/Figures/DataFig7.csv', header=True, index=None, sep='\t')
		fig8 = getFigDensityByLocation(df, path, 'Pseudogene')
		fig8.to_csv(path_or_buf=path+'/Results/Figures/DataFig8.csv', header=True, index=None, sep='\t')
		supfig1 = getFigBysubClass(df, path, 'Pseudogene')
		supfig1.to_csv(path_or_buf=path+'/Results/Figures/DataSupFig1.csv', header=True, index=None, sep='\t')

def build_arg_parser():
	parser = argparse.ArgumentParser(description = 'getDataFig')
	GITDIR = os.getcwd()+'/'
	parser.add_argument ('-p', '--path', default = GITDIR)
	return parser

if __name__ == '__main__':
	parser = build_arg_parser()
	arg = parser.parse_args()
	path = arg.path
	main(path)
