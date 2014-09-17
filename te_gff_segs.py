#!/usr/bin/env python

from optparse import OptionParser
import os, sys, glob, subprocess, tempfile
#import pdb

####################################################################
# te_gff_segs.py
#
# This script takes a GFF file folder as input and outputs a GFF file 
# in which all entries are broken into segments of 300 bp with an 
# overlap of 50 bp between segments. The breaking up of sequences is
# done on the basis of a MSA. Here're the steps involved:
#
# 1). For each repeat do the following:
#		-- Make a fasta file of all instances
#		-- Run a MSA on that FASTA file
#		-- Identify consensus columns in the MSA file
#		-- Divide the consensus into segments of 300 each with an 
#		   overlap of 50 bp
#		-- Divide instances into segments and output the GFF file 
#		   with segments
####################################################################
def main():
	usage = 'usage:%prog [options] <gff_folder> <msa_folder>'
	parser = OptionParser(usage)
	parser.add_option('-o', dest='output_dir', default='gff_segs_files')
	(options, args) = parser.parse_args()

	if len(args)!=2:
		parser.error('Must provide the GFF file and the folder with MSA')
	else:
		gff_folder = args[0]
		msa_folder = args[1]

	if os.path.isdir(options.output_dir):
		print >> sys.stderr, 'Writing output to %s' %(options.output_dir)
	else:
		os.mkdir(options.output_dir)
		print >> sys.stderr, 'Writing output to %s' %(options.output_dir)

	new_gff_files = glob.glob('/'.join([gff_folder, '*.gff']))
	for gff_file in new_gff_files:
		repeat = gff_file.split('.')[0].split('/')[-1]
		print repeat # Just to know where I am while executing the script
		msa_file = '/'.join([msa_folder, repeat + '.fa'])
		sequences = {}
		seg_gff = '/'.join([options.output_dir, repeat + '_seg.gff'])
		seg_gff_file = open(seg_gff, 'w')
		for line in open(msa_file,'r'):
			if line[0] == '>':
				header = line.strip()
			else:
				sequences[header] = line.strip()
		consensus_sequence = sequences['>Consensus']
		consensus_columns = []
		for i in range(0, len(consensus_sequence)):
			consensus_letter = consensus_sequence[i]
			if consensus_letter == 'x':
				consensus_columns.append(i)
		number_of_segments = 0 # Step 1: Find the number of segments in the consensus
		position = 0
		while position < len(consensus_columns):
			number_of_segments += 1
			position = (250*number_of_segments) + 50
		segments = {}
		i = 1
		segment_start = 0
		while i <= number_of_segments: # Divide the columns into segments
			segment = 'segment_' + str(i)
			segment_end = segment_start + 300
			segments[segment] = consensus_columns[segment_start:segment_end]
			segment_start = segment_end - 50
			i += 1
		segment_indices = {}
		for key in segments.keys(): # Get the start and end of each segment 
			segment_indices[key] = [segments[key][0], segments[key][-1]]
		for line in open(gff_file):
			contents = line.strip().split('\t')
			chromosome, start, end, strand, name = contents[0], int(contents[3]), int(contents[4]), contents[6], contents[8]
			family = name.split(';')[0]
			header = '>' + chromosome + ':' + str(start) + '-' + str(end)
			msa_sequence = sequences[header]
			msa_columns = {} # Map MSA columns to segments
			for i in range(0,len(msa_sequence)):
				if msa_sequence[i] not in ['.','-']:
					msa_columns[i] = []
			for column in msa_columns.keys():
				for key in segment_indices.keys():
					segment = key
					segment_start = segment_indices[key][0]
					segment_end = segment_indices[key][1]
					if column in range(segment_start, segment_end+1):
						if msa_columns[column] == []:
							msa_columns[column]=segment
						elif msa_columns[column]!= segment: # Since we're overlapping segments by 50 bp a column can be in more than 1 segment
							msa_columns[column]= msa_columns[column] + '+' + segment
			nucleotides = range(start,end+1)
			i = 0
			nucleotides_segments = {}
			feature_segments = []
			for key in sorted(msa_columns.keys()): # Map nucleotides to segments
				nucleotides_segments[nucleotides[i]] = msa_columns[key]
				column_segment = msa_columns[key].split('+')
				if len(column_segment) == 1:
					if column_segment[0] not in feature_segments:
						feature_segments.append(column_segment[0])
				else:
					if column_segment[0] not in feature_segments: 
						feature_segments.append(column_segment[0])
					if column_segment[1] not in feature_segments: # For the odd case when the first column of the MSA is in 2 segments
						feature_segments.append(column_segment[1])
				i+=1
			segments_nucleotides = {}
			for segment in feature_segments: # I need to initialize the hash table to empty lists for each segment so need feature_segments for keys.
				segments_nucleotides[segment] = []
			for key in sorted(nucleotides_segments.keys()): # Map segments to nucleotides
				nucleotide = key
				segment = nucleotides_segments[key].split('+')
				if len(segment)==1:
					segments_nucleotides[segment[0]].append(nucleotide)
				else:
					segment_1 = segment[0]
					segment_2 = segment[1]
					segments_nucleotides[segment_1].append(nucleotide)
					segments_nucleotides[segment_2].append(nucleotide)
			for key in segments_nucleotides.keys(): # Split the feature into several segments and print it to output!
				segment = key
				segment_start = segments_nucleotides[segment][0]
				segment_end = segments_nucleotides[segment][-1]
				segment_name = family + '; repeat "' + repeat + '_' + segment + '";' # Trying to recreate the naming pattern of repeat masker
				out_line = [chromosome, 'Dfam', 'Repeat', str(segment_start), str(segment_end), '.', strand, '.', segment_name]
				print >> seg_gff_file, '\t'.join(out_line)

	#seg_gff_files = '/'.join([options.output_dir, '*.gff'])
	#dfam_seg_gff_file = '/'.join([options.output_dir, 'dfam_seg.gff'])
	#subprocess.call('cat %s > %s' %(seg_gff_files, dfam_seg_gff_file), shell=True)
####################################################################
# main()
####################################################################
if __name__ == '__main__':
	main()
	#pdb.runcall(main)