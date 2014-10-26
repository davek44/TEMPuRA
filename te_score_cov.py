#!/usr/bin/env python

from optparse import OptionParser
import os, subprocess, sys, math, glob, numpy, pysam, stats, pdb, gff, time
import resource
#import pdb

#####################################################################
# te_score_cov.py
#
# This script inputs the BAM file and RM GFF file and outputs a BED 
# file with the position and score for each feature in the GFF file. 
# In this script, I compare and get the ratio of the number of reads 
# at the GFF file feature in the BAM file versus the control. The 
# log2 of this ratio is used as the score for that feature. Remember 
# to only include repeats with > 1000 occurrences and a certain 
# minimum coverage in the BED file.
#
# To Do:
#	- Reduce memory footprint (currently - 32 GB)
#####################################################################

def main():
	usage = 'usage:%prog [options] <bam_files> <gff_file>'
	parser = OptionParser(usage)
	parser.add_option('-o', dest='output_dir', default='te_score_cov', help='The directory in which you want to write the output files. [Default: %default]')
	parser.add_option('-c', dest='control_files', help='Control files for the input BAM files. [Default: %default]')
	parser.add_option('-l', dest='log', default=False, action='store_true', help='log2 coverage [Default: %default]')
	parser.add_option('-n', dest='min_num_repeats', default = -1, type='int', help='The minimum of instances of a repeat needed to report it as a seperate BED file. [Default: %default]')
	parser.add_option('-s', dest='min_score', type='float', default = -1.0, help='The minimum average score of a repeat needed to report it as a seperate BED file. [Default: %default]')
	(options, args) = parser.parse_args()

	if len(args)!=2:
		parser.error('Must provide both the bam file and gff file')
	else:
		bam_files = args[0]
		gff_file = args[1]

	bam_files = bam_files.split(',')
	if options.control_files:
		control_files = options.control_files.split(',')

	######################################################################
	# Check if output directory exists. If not, create it. Map the feature
	# ID's to the strand and name for each feature
	######################################################################

	if os.path.isdir(options.output_dir):
		print >> sys.stderr, 'Writing output to %s' %(options.output_dir)
	else:
		os.mkdir(options.output_dir)
		print >> sys.stderr, 'Writing output to %s' %(options.output_dir)

	gff_info = {}
	start_time = time.time()
	for line in open(gff_file,'r'):
		contents = line.strip().split('\t')
		chromosome, start, end = contents[0], int(contents[3]), int(contents[4])
		strand, name = contents[6], contents[8]
		feature_id = (chromosome, start, end)
		gff_info[feature_id] = (strand, name)

	time_diff = str(round((time.time()-start_time)/60,1))
	print 'Finshed mapping the feature ids to strand and name in %s minutes' % time_diff
	memory_in_use()

	######################################################################
	# Get the number of reads overlapping each feature in the gff file. 
	# Convert them to proper scores and print to a file.
	######################################################################

	###################################
	# Compute Coverage
	###################################

	start_time = time.time()
	coverage, events = compute_coverage(gff_file, bam_files)
	if options.control_files:
		coverage_control, events_control = compute_coverage(gff_file, control_files)

	time_diff = str(round((time.time()-start_time)/60,1))
	print >> sys.stderr, 'Finshed calculating coverage over all BAM files in %s minutes' % time_diff
	memory_in_use()

	###################################
	# Normalize
	###################################

	start_time = time.time()
	for feature_id in coverage:
		for i in range(len(coverage[feature_id])):
			coverage[feature_id][i] = (1+coverage[feature_id][i])/float(events)
			if options.control_files:
				coverage_control[feature_id][i] = (1+coverage_control[feature_id][i])/float(events_control)    

	time_diff = str(round((time.time()-start_time)/60,1))
	print >> sys.stderr, 'Finshed normalizing coverage over all feature ids in %s minutes' % time_diff
	memory_in_use()

	###################################
	# Get feature scores and write to
	# master score file
	###################################

	start_time = time.time()
	te_score_fd = '/'.join([options.output_dir,'te_score_file.txt'])
	te_score_file = open(te_score_fd, 'w')
	for feature_id in gff_info.keys():
		chromosome, start, end = feature_id
		strand, name = gff_info[feature_id]
		if feature_id not in coverage:
			score = 0
		else:
			feature_scores = []
			for i in range(len(coverage[feature_id])):
				if options.log:
					cov = math.log(coverage[feature_id][i], 2)
				else:
					cov = coverage[feature_id][i]
				if options.control_files:
					if options.log:
						cov -= math.log(coverage_control[feature_id][i], 2)
					else:
						cov = cov / coverage_control[feature_id][i]
				feature_scores.append(cov)
			score = numpy.mean(feature_scores)
		out_line = '\t'.join([chromosome, str(start), str(end), name, str(score), strand])
		print >> te_score_file, out_line

	te_score_file.close()
	time_diff = str(round((time.time()-start_time)/60,1))
	print >> sys.stderr, 'Finshed calculating scores for all feature ids in %s minutes' % time_diff
	print >> sys.stderr, 'Master score file written to %s' %(te_score_fd)
	memory_in_use()

	######################################################################
	# Filter the score BED file and only keep the repats with a minimum
	# number of instances and an average minimum score.
	######################################################################

	least_score = options.min_score
	least_lines = options.min_num_repeats
	
	te_score_dir = options.output_dir
	split_master(te_score_fd, te_score_dir, 'repeats', least_score, least_lines)
	memory_in_use()

	######################################################################
	# Filter the score BED file and only keep the repats with a minimum
	# number of instances and an average minimum score.
	######################################################################

	split_master(te_score_fd, te_score_dir, 'families', least_score, least_lines)
	print >> sys.stderr, 'Script finished. Check folder - %s for repeats BED files' %(options.output_dir)

################################################################################
# compute_coverage
#
# Input
#  gff_file:    GFF file of genome features.
#  bam_files:   BAM or GFF files of reads alignments.
#
# Output
#  coverage:      Dict mapping feature_id's to coverage arrays.
#  events:        Total number of events.
################################################################################

def compute_coverage(gff_file, bam_files):
	coverage = {}
	for line in open(gff_file):
		a = line.split('\t')
		chrom = a[0]
		start = int(a[3])
		end = int(a[4])
		feature_length = abs(end-start+1)
		feature_id = (chrom, start, end)
		if not feature_id in coverage:
			coverage[feature_id] = [0]*feature_length

	events = 0
	for bam_file in bam_files:
		print >> sys.stderr, 'Computing coverage for %s' % bam_file
		multi_maps = {}
		for aligned_read in pysam.Samfile(bam_file, 'rb'):
			try:
				nh_tag = aligned_read.opt('NH')
			except:
				nh_tag = 1.0

			if aligned_read.is_paired:
				events += 0.5/nh_tag
			else:
				events += 1.0/nh_tag

			if nh_tag > 1:	
				multi_maps[aligned_read.qname] = nh_tag

		p = subprocess.Popen('intersectBed -sorted -wo -bed -abam %s -b %s' % (bam_file, gff_file), shell=True, stdout=subprocess.PIPE)
		for line in p.stdout:
			a = line.split('\t')
			rstart = int(a[1])+1  # convert back to 1-based gff from bed
			rend = int(a[2])
			rheader = a[3]
			if rstart < rend:
				acol = 12
				achrom = a[acol]
				astart = int(a[acol+3])
				aend = int(a[acol+4])
				astrand = a[acol+6]
				feature_id = (achrom, astart, aend)
				cov_start = max(rstart, astart)
				cov_end = min(rend, aend)
				alength = aend - astart + 1
				if astrand == '+':
					inc_start = (cov_start - astart) 
					inc_end = (cov_end - astart + 1) 
				else:
					inc_start = (aend - cov_end) 
					inc_end = (aend - cov_start + 1) 
				if inc_start != None:
					if rheader in multi_maps:
						mm = multi_maps[rheader]
					else:
						rheader_base = rheader[:rheader.rfind('/')]
						if rheader_base in multi_maps:
							mm = multi_maps[rheader_base]
						else:
							mm = 1.0
					for i in range(inc_start, inc_end):
						coverage[feature_id][i] += 1.0/mm
		p.communicate()

	return coverage, events

################################################################################
# split_master
#
# Input
#  master_file_name:	A master file with scores for each repeat
#  identifier:			A string to decide whether to split families/repeats
#  least_score:			Float or int deciding threshold mean score
#  least_lines:			Int deciding threshold number of lines 
#
# This function splits the master file based on the identifier into 1 file for 
# each indentifier and prints a summary file. I also remove the individual BED
# files which fail to meet the criteria of either least_score or least_lines.
# If those parameters are not provided, based on the mean I guess them.
################################################################################

def split_master(master_file_name, master_dir, identifier, least_score, least_lines):

	if identifier == 'repeats':
		split_contents = 1
		summary_file = 'summary_statistics_repeats.txt'
		sorted_summary_file = 'sorted_summary_statistics_repeats.txt'
	elif identifier == 'families':
		split_contents = 0
		summary_file = 'summary_statistics_families.txt'
		sorted_summary_file = 'sorted_summary_statistics_families.txt'
	else:
		exit(1)

	output_dir = '/'.join([master_dir, identifier])
	if os.path.isdir(output_dir):
		print >> sys.stderr, 'Writing family families to %s' %(output_dir)
	else:
		os.mkdir(output_dir)
		print >> sys.stderr, 'Writing output to %s' %(output_dir)

	###################################
	# Split the score's file into 1 
	# file for every repeat
	###################################

	repeat_families = {}
	for line in open(master_file_name, 'r'):
		contents = line.split('\t')
		repeat_family = contents[3].split(';')[split_contents].split('"')[1]
		if len(repeat_family.split('/')) > 1:
			repeat_family = repeat_family.split('/')[1]
		if repeat_family not in repeat_families:
			repeat_family_file = '/'.join([output_dir, repeat_family]) + '.bed'
			repeat_families[repeat_family] = open(repeat_family_file, 'w')
		
		print >> repeat_families[repeat_family], line.strip()
		
	#print >> sys.stderr, 'Split the master score file into 1 BED file for each repeat family'

	for key in repeat_families.keys():
		repeat_families[key].close()
	
	###################################
	# Write the Summary Statistics
	# file
	###################################

	bed_files = '/'.join([output_dir, '*.bed'])
	bed_file_names = glob.glob(bed_files)
	summary_statistics = '/'.join([output_dir, summary_file])
	summary_statistics_file = open(summary_statistics, 'w')
	
	print >> summary_statistics_file, identifier + '\tMean\tMedian\tMinimum\tMaximum\tStdev\tVariance\tLines'
	
	all_families_score = []
	all_families_lines = []
	for bed_file_name in bed_file_names:
		repeat = bed_file_name.split('.')[0].split('/')[-1]
		scores = []
		num_lines = 0
		bed_file = open(bed_file_name, 'r')
		for line in bed_file:
			score = float(line.split('\t')[4])
			scores.append(score)
			num_lines +=1
		bed_file.close()
		mean_score = numpy.mean(scores)
		all_families_score.append(mean_score)
		all_families_lines.append(num_lines)
		score_variance = sum((mean_score - value)**2 for value in scores)/len(scores)
		min_score = min(scores)
		max_score = max(scores)
		score_stdev = numpy.std(scores)
		median_score = numpy.median(scores)
		out_line = '\t'.join([repeat, str(mean_score), str(median_score), str(min_score), str(max_score), str(score_stdev), str(score_variance), str(num_lines)])
		print >> summary_statistics_file, out_line
		
	summary_statistics_file.close()
	print >> sys.stderr, 'Finished making the %s Summary Statistics file. File written to %s' %(identifier, summary_statistics)

	###################################
	# If -n and -s are not provided try
	# to guess appropriate values
	###################################

	mean_families_score = numpy.mean(all_families_score)
	mean_families_lines = numpy.mean(all_families_lines)
	if least_score == -1:
		min_score = mean_families_score
	else:
		min_score = least_score
	if least_lines == -1:
		min_num_repeats = mean_families_lines
	else:
		min_num_repeats = least_lines

	###################################
	# Remove the BED files not matching
	# the minimum mean score and number
	# of lines
	###################################

	sorted_summary_statistics = '/'.join([output_dir, sorted_summary_file])
	subprocess.call("(head -n 1 %s && tail -n +2 %s | sort -t$'\t' -nr -k2,6) > %s" %(summary_statistics, summary_statistics, sorted_summary_statistics), shell=True)
	for line in open(sorted_summary_statistics,'r'):
		if line.strip()[-5:] != 'Lines':
			contents = line.strip().split('\t')
			num_lines = int(contents[-1])
			mean_score = float(contents[1])
			repeat_family = contents[0]
			bed_file_name = '/'.join([output_dir, repeat_family]) + '.bed'
			if mean_score < min_score:
				os.remove(bed_file_name)
			elif num_lines < min_num_repeats:
				os.remove(bed_file_name)

	print >> sys.stderr, 'Removed the %s with average score less than %f' %(identifier, round(min_score,2))
	print >> sys.stderr, 'Removed the %s with less than %d instances' %(identifier, round(min_num_repeats,2))
################################################################################
# memory_in_use
# 
# Input:	No input arguments!
# Output:	The total memory usage for calling process (Resident Set Size - RSS)
# 
# A small function to print the total memory being used by the process in GB
################################################################################

def memory_in_use():
	kb_memory_used = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	gb_memory_used = kb_memory_used/(1024*1024)
	print >> sys.stderr, 'Current memory being used is %d kb or %f GB' % (kb_memory_used, gb_memory_used)

#####################################################################
# main()
#####################################################################
if __name__ == '__main__':
 	main()
 	#pdb.runcall(main)