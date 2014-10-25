#!/usr/bin/env python

from optparse import OptionParser
import subprocess, tempfile

################################################################################
# te_cov.py
#
# plot occurrences of a known motif along a repeat, so we can compare to the 
# regression coefficients and mutual information
################################################################################

################################################################################
# main
################################################################################
def main():
	usage='usage:%prog [options] <msa_file>'
	parser = OptionParser(usage)
	parser.add_option('-l', dest='possum_library', help='PWM library to use for matching motifs')
	(options, args) = parser.parse_args()

	if len(args)!=1:
		parser.error('Must provide MSA FASTA File. %s' % usage)
	else:
		msa_file = args[0]

    ############################################
    # Read all the MSA sequences and print them
    # to a temporary file after removing gaps.
    ############################################
    msa_sequences = {}
    for line in open(msa_file):
    	if line[0] == '>':
    		header = line.strip()[1:]
    		msa_sequences[header] = ''
    	else:
    		msa_sequences[header] += line.strip()

    msa_sequences.pop('Consensus') # Remove Consensus sequence
    msa_sequences_fd, msa_sequences_file_name = tempfile.mkstemp()
	msa_sequences_file = open(msa_sequences_file_name, 'w')
	for header in msa_sequences.keys():
		print >> msa_sequences_file, header
		no_gaps_sequence = msa_sequences[header].translate(None, '.-') # Remove Gaps
		print >> msa_sequences_file, no_gaps_sequence

	msa_sequences_file.close()

	############################################
	# Create a PoSSuM Index from the temporary
	# file created above
	############################################
	subprocess.call('mkvtree -db %s -indexname msa_index -dna -tis -suf -lcp -skp -v', shell=True)
	subprocess.call('possumfreqs -db %s > frequencies.txt' % ,shell=True)
	subprocess.call('possumdist -pr %s -dna -freq frequencies.txt -pdis dist.gz' % options.possum_library, shell=True)
	subprocess.call('possumsearch -pr %s -dna -db %s -freq frequencies.txt -lazy -esa -pval 1e-6 -fn -rc -format tabs > msa_motifs.txt' %(options.possum_library, msa_sequences_file_name), shell=True)
	subprocess.call('possum2gff.py msa_motifs.txt > msa_motifs.gff', shell=True)

################################################################################
# __main__
################################################################################
if __name__ == '__main__':
	main()


