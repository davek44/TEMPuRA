#!/usr/bin/env python

from optparse import OptionParser
import ggplot

##################################################################################
# te_mut_info.py
#
# This script inputs a MSA file and a BED file with scores. The output is mutual
# information between each column of the consensus in MSA file and the scores.
# 
# To Do:
#     - Currently, I discretize my scores before computing MI. It would be good to 
#    calculate MI for continuous scores.
##################################################################################
def main():
    usage='usage:%prog [options] <bed_file> <msa_file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='consensus_pct', default=0.5, type='float', help='Required proportion of columns with a valid nt to consider it a consensus column [Default: %default]')
    parser.add_option('-d', dest='dfam_consensus', action='store_true', help='Pass the option if you want to use Consensus as defined by Dfam')
    #parser.add_option('-j', dest='condense_pct', type='float', help='Required proportion of entries to be same between 2 columns for them to be merged')
    #parser.add_option('-n', dest='discretize_bins', type='int', help='The number of bins you want to discretize the scores into')
    parser.add_option('-o', dest='output_pre', type='string', help='Prefix of the output files')
    (options, args) = parser.parse_args()

    if len(args)!=2:
        parser.error('Must provide both the BED file and MSA file. Check %s' %usage)
    else:
        bed_file = args[0]
        msa_fasta_file = args[1]

    ##################################################
    # hash scores
    ##################################################
    seq_scores = {}
    for line in open(bed_file):
        a = line.split('\t')
        header = a[0] + ':' + a[1] + '-' + a[2]
        score = float(a[4])
        seq_scores[header] = score

    ##################################################
    # define consensus
    # define columns to condense for regression
    ##################################################
    msa_sequences = {}
    for line in open(msa_fasta_file):
        if line[0] == '>':
            header = line.strip()
            msa_sequences[header] = ''
        else:
            msa_sequences[header] += line.strip()

    if options.dfam_consensus is True:
        consensus_sequence = msa_sequences.pop('>Consensus')
        sequence_length = len(consensus_sequence)
        consensus_columns = []
        for i in range(0,len(consensus_sequence)):
            if consensus_sequence[i] == 'x':
                consensus_columns.append(i)
    else:
        consensus_columns = define_consensus(msa_fasta_file, options.consensus_pct)
        #sample_sequence = msa_sequences.pop('>Consensus')
        #sequence_length = len(sample_sequence)

    #hamming_cutoff = int(sequence_length - options.condense_pct*sequence_length)
    #condensed_columns, columns_ls_remove = column_condense(msa_sequences, consensus_columns, hamming_cutoff)

    ##################################################
    # map sequences to feature vectors
    ##################################################
    # initialize the dictionary with score and position/nt features
    df_mi = {'Score':[]}
    for i in range(len(consensus_columns)):
        position = i+1
        df_mi[position] = []

    header = ''
    for line in open(msa_fasta_file):
        if line[0] == '>':
            if header and header != 'Consensus':
                # process seq
                df_mi['Score'].append(seq_scores[header])
                for i in range(len(consensus_columns)):
                    position = i+1
                    seq_i = consensus_columns[i]
                    nt = seq[seq_i].upper()
                    df_mi[position].append(nt)
            
            header = line[1:].rstrip()
            seq = ''

        else:
            seq += line.rstrip()
    
    if header and header != 'Consensus':
        # process last seq
        df_mi['Score'].append(seq_scores[header])
        for i in range(len(consensus_columns)):
            position = i+1
            seq_i = consensus_columns[i]
            nt = seq[seq_i].upper()
            df_mi[position].append(nt)

    r_script = '/Users/chinmayshukla/Documents/Research/Sandbox-TEMPuRA/bin/r/te_mut_info.r'
    ggplot.plot(r_script, df_mi, [options.output_pre])
################################################################################
# define_consensus
#
# Input
#  msa_fasta_file:
#  consensus_pct:    Float above which we consider the column to be consensus.
#
# Output
#   consensus_cols:  List of consensus column indexes.
################################################################################
def define_consensus(msa_fasta_file, consensus_pct):
    valid_nts = ['A','C','G','T']
    
    column_counts = []
    header = ''
    seq_count = 0
    for line in open(msa_fasta_file):
        if line[0] == '>':            
            if header and header != 'Consensus':  # avoid DFAM
                seq_count += 1
                for i in range(len(seq)):
                    if seq[i].upper() in valid_nts:
                        while i >= len(column_counts):
                            column_counts.append(0)
                        column_counts[i] += 1

            header = line[1:].rstrip()
            seq = ''
            
        else:
            seq += line.rstrip()
                        
    if header and header != 'Consensus':
        seq_count += 1
        for i in range(len(seq)):
            if seq[i].upper() in valid_nts:
                while i >= len(column_counts):
                    column_counts.append(0)
                column_counts[i] += 1

    consensus_columns = []
    for i in range(len(column_counts)):
        if column_counts[i] / float(seq_count) > consensus_pct:
            consensus_columns.append(i)

    return consensus_columns

##################################################################################
# main()
##################################################################################
if __name__ == '__main__':
    main()