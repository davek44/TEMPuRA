#!/usr/bin/env python

import random, sys

################################################################################
# generate_seq_scores.py
#
# Generates random sequences, with a functional motif at the start. Randomly
# score sequences from a pair of gaussian distributions, according to whether
# the functional motif was altered.
################################################################################
def main():
    
    bed_out = open('test_scores.bed','w')
    msa_out = open('test_msa.fa','w')

    #random.seed(1234)
    nucleotides = ['A','C','T','G']

    num_seqs = 2000
    seq_length = 10
    motif_length = 5
    mut_rate = 0.3

    # randomly choose a consensus sequence
    consensus_sequence = []
    for i in range(0,seq_length):
        nucleotide = nucleotides[random.randint(0,3)]
        consensus_sequence.append(nucleotide)
    consensus_motif = ''.join(consensus_sequence[:motif_length])
    print consensus_motif

    for n in range(num_seqs):
        # mutate sequence
        seq = ''
        for i in range(len(consensus_sequence)):
            if random.random() < mut_rate:
                seq += nucleotides[random.randint(0,3)]
            else:
                seq += consensus_sequence[i]

        # print FASTA
        header = 'seq%d' % (n+1)
        print >> msa_out, '>%s\n%s' % (header, seq)

        # check motif
        seq_motif = seq[:motif_length]
        if seq_motif == consensus_motif:
            raw_score = 10
        else:
            raw_score = 0

        # randomize score
        score = raw_score + random.gauss(0,1)

        # print BED
        fake_start = n*seq_length
        fake_end = fake_start + seq_length
        cols = ['chrX', str(fake_start), str(fake_end), header, str(score), '+']
        print >> bed_out, '\t'.join(cols)

    bed_out.close()
    msa_out.close()

    
################################################################################
# main()
################################################################################
if __name__ == '__main__':
	main()
