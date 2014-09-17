#!/usr/bin/env python

from optparse import OptionParser
import os, sys, subprocess, tempfile, glob, numpy
#import pdb

#####################################################################
# te_score_peaks.py
#
# This script inputs the peaks BED file and RM GFF file and outputs a
# BED file with the position and score for each feature in the GFF 
# file. In this script, I use peaks BED file made using AREM and check
# if the feature in the RM GFF file intersects with the peaks BED file;
# hence assigning each feature a score of 0 or 1. Remember to only 
# include repeats with > 1000 occurrences and a certain minimum number 
# of peaks in the BED file.
#####################################################################
def main():
	usage='usage:%prog [options] <peaks_file> <gff_file>'
	parser = OptionParser(usage)
	parser.add_option('-o', dest='output_dir', default='te_score_peaks', help='The directory in which you want to write the output files. [Default: %default]')
	parser.add_option('-n', dest='num_repeats', type='int', default=1, help='The minimum of instances of a repeat needed to report it as a seperate BED file. [Default: %default]')
	parser.add_option('-s', dest='min_score', type='float', default=0.0, help='The minimum average score of a repeat needed to report it as a seperate BED file. [Default: %default]')
	(options, args) = parser.parse_args()

	if len(args) != 2:
		parser.error('Must provide the peaks BED file and GFF file')
	else:
		peaks_file = args[0]
		gff_file = args[1]

	###################################
	# Step 0: Check if output directory
	# exists. If not, create it.
	###################################
	
	if os.path.isdir(options.output_dir):
		print >> sys.stderr, 'Writing output to %s' %(options.output_dir)
	else:
		os.mkdir(options.output_dir)
		print >> sys.stderr, 'Writing output to %s' %(options.output_dir)

	peaks_features_file = '/'.join([options.output_dir, 'peaks_features.txt'])
	subprocess.call('intersectBed -a %s -b %s -wo > %s' %(peaks_file, gff_file, peaks_features_file), shell=True)
	print >> sys.stderr, 'Finished intersecting the peaks file and gff file'

	###################################
	# Step 1: For each entry in the GFF
	# file assign a score (0 or 1) and 
	# write to a new scores master file
	###################################

	peaks_features = {}
	for line in open(peaks_features_file):
		contents = line.strip().split('\t')
		chromosome, start, end, strand = contents[5], contents[8], contents[9], contents[11]
		feature_id = '\t'.join([chromosome, start, end, strand])
		peaks_features[feature_id] = 1

	te_score_fd = '/'.join([options.output_dir,'te_score_file.txt'])
	te_score_file = open(te_score_fd, 'w')
	for line in open(gff_file):
		contents = line.strip().split('\t')
		chromosome, start, end, strand, name = contents[0], contents[3], contents[4], contents[6], contents[8]
		feature_id = '\t'.join([chromosome, start, end, strand])
		if feature_id in peaks_features:
			score = 1
		else:
			score = 0
		out_line = '\t'.join([chromosome, start, end, strand, name, str(score)])
		print >> te_score_file, out_line
	te_score_file.close()
	print >> sys.stderr, 'Master score file written to %s' %(te_score_fd)

	###################################
	# Step 2: Filter the score BED file
	# and only keep the repats with a 
	# minimum number of instances and an
	# average minimum score.
	###################################

	repeats = {}
	for line in open(te_score_fd, 'r'):
		contents = line.split('\t')
		repeat = contents[4].split(';')[1].split('"')[1]
		if repeat not in repeats:
			repeat_file = '/'.join([options.output_dir, repeat]) + '.bed'
			repeats[repeat] = open(repeat_file, 'w')
		
		print >> repeats[repeat], line.strip()
		
	print >> sys.stderr, 'Split the master score file into 1 BED file for each repeat'

	for key in repeats.keys():
		repeats[key].close()

	bed_files = '/'.join([options.output_dir, '*.bed'])
	bed_files = glob.glob(bed_files)
	summary_statistics = '/'.join([options.output_dir, 'summary_statistics.txt'])
	summary_statistics_file = open(summary_statistics, 'w')
	print >> summary_statistics_file, 'Repeat\tMean\tMedian\tMinimum\tMaximum\tStdev\tVariance\tLines'
	for bed_file in bed_files:
		repeat = bed_file.split('.')[0].split('/')[-1]
		scores = []
		num_lines = 0
		for line in open(bed_file, 'r'):
			score = float(line.split('\t')[-1])
			scores.append(score)
			num_lines +=1
		mean_score = numpy.mean(scores)
		score_variance = sum((mean_score - value)**2 for value in scores)/len(scores)
		min_score = min(scores)
		max_score = max(scores)
		score_stdev = numpy.std(scores)
		median_score = numpy.median(scores)
		out_line = '\t'.join([repeat, str(mean_score), str(median_score), str(min_score), str(max_score), str(score_stdev), str(score_variance), str(num_lines)])
		print >> summary_statistics_file, out_line
		if mean_score < options.min_score:
			subprocess.call('rm %s' %(bed_file), shell=True)

	print >> sys.stderr, 'Wrote summary statistics to %s' %(summary_statistics)
	print >> sys.stderr, 'Removed the repeats with average score less than %f' %(round(options.min_score,2))

	bed_files = '/'.join([options.output_dir, '*.bed'])
	bed_files = glob.glob(bed_files)
	for bed_file in bed_files:
		p = subprocess.Popen('wc -l %s' %(bed_file), shell=True, stdout=subprocess.PIPE)
		for line in p.stdout:
			num_repeats = int(line.strip().split(' ')[0])
		if num_repeats < options.num_repeats:
			subprocess.call('rm %s' %(bed_file), shell=True)

	print >> sys.stderr, 'Removed the repeats with less than %d instances' %(round(options.num_repeats,2))
	print >> sys.stderr, 'Script finished. Check folder - %s for repeats BED files' %(options.output_dir)
#####################################################################
# main()
#####################################################################
if __name__ == '__main__':
	main()
	#pdb_runcall(main)