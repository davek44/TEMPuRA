#!/usr/bin/env python

from optparse import OptionParser
import pandas as pd
import statsmodels.api as sm
import os, sys, subprocess, ggplot, collections

#####################################################################
# te_feature_weights.py
#
# This script inputs a bed file describing various features and their
# scores. The output is a file with weights for each nucleotide in 
# the input bed file sequence.
#
# For a gap, should I give weights for 3 nucleotides or 4 nucleotides?
#####################################################################
def main():
	usage = 'usage:%prog [options] <bed_file> <msa_file>'
	parser = OptionParser(usage)
	parser.add_option('-m', dest='model_output_file', default='ols_summary.txt', help='The file to write the model summary')
	parser.add_option('-p', dest='plot_output_file', default='weights_plot.pdf', help='The file to print the plot of index versus weight')
	(options, args) = parser.parse_args()

	if len(args)!=2:
		parser.error('Must provide bed_file with scores')
	else:
		bed_file = args[0]
		msa_file = args[1]

	bed_features_sequences = {}
	for line in open(msa_file):
		if line[0]=='>':
			header = line.strip()
			bed_features_sequences[header] = ''
		else:
			bed_features_sequences[header] += line.strip()

	consensus_sequence = bed_features_sequences['>Consensus']
	consensus_columns = []
	for i in range(0,len(consensus_sequence)):
		if consensus_sequence[i] == 'x':
			consensus_columns.append(i)
	quaternary_conversion_dict = {'A':[1,0,0], 'a':[1,0,0], 'C':[0,1,0], 'c':[0,1,0], 'T':[0,0,0], 't':[0,0,0], 'G':[0,0,1], 'g':[0,0,1], '.':[0.25,0.25,0.25], '-':[0.25,0.25,0.25]}
	df_dict = {}
	df_dict['Score'] = []
	for i in range(0, len(consensus_columns)):
		position = str(i+1)
		df_dict[position+'_A'] = []
		df_dict[position+'_C'] = []
		#df_dict[position+'_T'] = []
		df_dict[position+'_G'] = []

	for line in open(bed_file):
		contents = line.strip().split('\t')
		score = float(contents[-1])
		df_dict['Score'].append(score)
		chromosome, start, end, strand = contents[0], int(contents[1]), int(contents[2]), contents[3] 
		for key in bed_features_sequences.keys():
			header = key
			if header != '>Consensus':
				header_chromosome = header.split(':')[0][1:]
				header_position = header.split(':')[1].split('-')
				header_start, header_end = int(header_position[0]), int(header_position[1])
				header_range = range(header_start, header_end+1)
				if (header_chromosome == chromosome and start in header_range and end in header_range):
					feature_sequence = bed_features_sequences[header]
		j = 1
		for i in range(0, len(feature_sequence)):
			if i in consensus_columns:
				position = str(j)
				letter = feature_sequence[i]
				letter_score = quaternary_conversion_dict[letter]
				df_dict[position+'_A'].append(letter_score[0])
				df_dict[position+'_C'].append(letter_score[1])
				#df_dict[position+'_T'].append(letter_score[2])
				df_dict[position+'_G'].append(letter_score[2])
				j += 1

	df = pd.DataFrame(df_dict)
	score = df['Score']
	X = df.drop('Score', axis=1)
	print >> sys.stderr, 'Read in all the sequences and scores. Now fitting the model'
	mod = sm.OLS(score, X)
	res = mod.fit()
	model_output_file = open(options.model_output_file,'w')
	print >> model_output_file, res.summary()
	model_output_file.close()
	print >> sys.stderr, 'Fit an OLS model and print the summary to %s' %(options.model_output_file)

	position_weights = collections.defaultdict(list)
	#for i in range(1,j):
	#	position_weights[i] = []
	flag = False
	for line in open(options.model_output_file, 'r'):
		if line[0:2] == '==':
			flag = False
		elif line[0:2] == '--':
			flag = True
		elif flag:
			contents = line.split()
			position = int(contents[0].split('_')[0])
			weight = float(contents[1])
			position_weights[position].append(weight)

	df_dict = {'Position':[], 'Nucleotide':[], 'Weight':[]}
	print '\t'.join(df_dict.keys())
	for position in position_weights.keys():
		weight_A, weight_C, weight_G = position_weights[position]
		weight_T = 0.0
		nucleotide_weights = [weight_A, weight_C, weight_T, weight_G]
		nucleotide_order = ['A','C','T','G']
		min_weight = min(nucleotide_weights)
		for i in range(0, len(nucleotide_weights)):
			nucleotide_weights[i] = nucleotide_weights[i] - min_weight
		max_weight = max(nucleotide_weights)
		for i in range(0, len(nucleotide_weights)):
			nucleotide_weights[i] = nucleotide_weights[i]/max_weight
		for i in range(0, len(nucleotide_weights)):
			df_dict['Position'].append(position)
			df_dict['Nucleotide'].append(nucleotide_order[i])
			df_dict['Weight'].append(nucleotide_weights[i])
			print '\t'.join([str(position), nucleotide_order[i], str(nucleotide_weights[i])])

	print >> sys.stderr, 'Now plotting the weights of different nucleotides along each position'
	ggplot.plot('/Volumes/Scratch/Chinmay/TEMPuRA/plot_te_weights_heatmap.r', df_dict, [options.plot_output_file])
	print >> sys.stderr, 'All Done. Check output files'

#####################################################################
# main()
#####################################################################
if __name__ == '__main__':
	main()