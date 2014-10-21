#!/usr/bin/env python
from optparse import OptionParser
import pandas as pd
import statsmodels.api as sm
import os, sys, subprocess, collections
import ggplot
import tempura

################################################################################
# te_score_plots.py
#
# Given a BED file containing scores and a MSA fasta file, perform a regression
# on the nucleotides at each position in the MSA against the scores.
################################################################################


################################################################################
# main
################################################################################
def main():
    usage = 'usage:%prog [options] <bed_file> <msa_file>'
    parser = OptionParser(usage)
    parser.add_option('-c', dest='consensus_pct', default=0.5, type='float', help='Required proportion of columns with a valid nt to consider it a consensus column [Default: %default]')
    parser.add_option('-m', dest='model_output_file', default='ols_summary.txt', help='The file to write the model summary')
    parser.add_option('-p', dest='plot_output_file', default='weights_plot.pdf', help='The file to print the plot of index versus weight')
    (options, args) = parser.parse_args()
    
    if len(args) != 2:
        parser.error('Must provide BED file with scores and MSA fasta file')
    else:
        bed_file = args[0]
        msa_fasta_file = args[1]

    ##################################################
    # hash scores
    ##################################################
    seq_scores = {}
    for line in open(bed_file):
        a = line.split('\t')
        header = a[3]
        score = float(a[4])
        seq_scores[header] = score

        
    ##################################################
    # define consensus
    ##################################################
    consensus_columns = define_consensus(msa_fasta_file, options.consensus_pct)


    ##################################################
    # map sequences to feature vectors
    ##################################################
    quaternary_conversion_dict = {'A':[1,0,0], 'C':[0,1,0], 'G':[0,0,1], 'T':[0,0,0], 'N':[0.25,0.25,0.25], '.':[0.25,0.25,0.25], '-':[0.25,0.25,0.25]}

    # initialize the dictionary with score and position/nt features
    df_dict = {'Score':[]}
    for i in range(len(consensus_columns)):
        position = str(i+1)
        df_dict[position+'_A'] = []
        df_dict[position+'_C'] = []
        df_dict[position+'_G'] = []

    header = ''
    for line in open(msa_fasta_file):
        if line[0] == '>':
            if header and header != 'Consensus':
                # process seq
                df_dict['Score'].append(seq_scores[header])
                for i in range(len(consensus_columns)):
                    position = str(i+1)
                    seq_i = consensus_columns[i]
                    nt = seq[seq_i].upper()
                    nt_conv = quaternary_conversion_dict[nt]
                    df_dict[position+'_A'].append(nt_conv[0])
                    df_dict[position+'_C'].append(nt_conv[1])
                    df_dict[position+'_G'].append(nt_conv[2])
            
            header = line[1:].rstrip()
            seq = ''

        else:
            seq += line.rstrip()
    
    if header and header != 'Consensus':
        # process last seq
        df_dict['Score'].append(seq_scores[header])
        for i in range(len(consensus_columns)):
            position = str(i+1)
            seq_i = consensus_columns[i]
            nt = seq[seq_i].upper()
            nt_conv = quaternary_conversion_dict[nt]
            df_dict[position+'_A'].append(nt_conv[0])
            df_dict[position+'_C'].append(nt_conv[1])
            df_dict[position+'_G'].append(nt_conv[2])


    ##################################################
    # perform learning
    ##################################################
    # add y-intercept term
    df_dict['Const'] = [1]*len(df_dict['Score'])

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


    ##################################################
    # read output
    ##################################################
    position_weights = collections.defaultdict(list)
    flag = False
    for line in open(options.model_output_file, 'r'):
        if line[0:2] == '==':
            flag = False
        elif line[0:2] == '--':
            flag = True
        elif flag:
            contents = line.split()
            if contents[0] != 'Const':
                position = int(contents[0].split('_')[0])
                weight = float(contents[1])
                position_weights[position].append(weight)

    df_dict = {'Position':[], 'Nucleotide':[], 'Weight':[]}
    #print '\t'.join(df_dict.keys())
    for position in position_weights.keys():
        weight_A, weight_C, weight_G = position_weights[position]
        weight_T = 0.0
        nucleotide_weights = [weight_A, weight_C, weight_T, weight_G]
        nucleotide_order = ['A','C','T','G']

        min_weight = min(nucleotide_weights)
        for i in range(0, len(nucleotide_weights)):
            nucleotide_weights[i] = nucleotide_weights[i] - min_weight

        for i in range(0, len(nucleotide_weights)):
            df_dict['Position'].append(position)
            df_dict['Nucleotide'].append(nucleotide_order[i])
            df_dict['Weight'].append(nucleotide_weights[i])
            #print '\t'.join([str(position), nucleotide_order[i], str(nucleotide_weights[i])])

    print >> sys.stderr, 'Now plotting the weights of different nucleotides along each position'
    #ggplot.plot('%s/te_score_plots.r' % tempura.r_dir, df_dict, [options.plot_output_file])
    ggplot.plot('/Users/chinmayshukla/Documents/Research/TEMPuRA/r/te_score_plots.r', df_dict, [options.plot_output_file])
    print >> sys.stderr, 'All Done. Check output files'


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


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
