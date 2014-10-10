#!/usr/bin/env python
from optparse import OptionParser
import copy, glob, os, random, re, shutil, stats, subprocess, sys, pdb, tempfile
import pysam
import dna, gff, ggplot, te, util
import tempura

################################################################################
# te_cov.py
#
# Plot BAM coverage across a TE family.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <dfam_te> <bam1> ... <bam2> ... <bamN>'
    parser = OptionParser(usage)

    parser.add_option('-l', dest='labels', help='BAM labels (comma separated) for plots')
    parser.add_option('-s', dest='stranded', default=False, action='store_true', help='Seq is stranded so consider orientation [Default: %default]')
    parser.add_option('-t', dest='threads', default=1, type='int', help='Number of threads to run for multi-TE mode [Default: %default]')

    # MSA
    parser.add_option('-c', dest='clean', default=False, action='store_true', help='Clean up removing temp files like the MSA computed for the reads [Default: %default]')
    parser.add_option('-m', dest='msa_saved', default=False, action='store_true', help='MSAs computed by HMMer are already done and in the output directory [Default: %default]')

    # limitations
    parser.add_option('-i', dest='min_cov', default=100, type='int', help='Minimum number of reads/nt to proceed after intersectBed [Default: %default]')
    parser.add_option('-a', dest='max_reads', default=10000, type='int', help='Maximum number of reads to hmmalign [Default: %default]')

    parser.add_option('-o', dest='out_dir', default='te_cov', help='Output directory [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) < 2:
        parser.error('Must provide DFAM TE family and BAM')
    else:
        dfam_te = args[0]
        bam_files = args[1:]

    ############################################
    # launch recursively
    ############################################
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    if dfam_te == '.':
        processes = []
        for te_gff in glob.glob('%s/gffs/*.gff' % tempura.dfam_dir):
            te = os.path.splitext(os.path.split(te_gff)[1])[0]

            cmd_argv = copy.copy(sys.argv)

            # change out_dir
            i = 0
            while i < len(cmd_argv) and cmd_argv[i] != '-o':
                i += 1

            te_out = '%s/%s' % (options.out_dir,te)
            if i == len(cmd_argv):
                cmd_argv = cmd_argv[0] + ['-o', te_out] + cmd_argv[1:]
            else:
                cmd_argv[i+1] = te_out

            # change dfam_te
            cmd_argv[-len(bam_files)-1] = te

            cmd = ' '.join(cmd_argv)
            processes.append(cmd)

        util.exec_par(processes, options.threads, print_cmd=True)

        exit()

    ############################################
    # prep
    ############################################
    if options.labels:
        labels = options.labels.split(',')
    else:
        labels = []

    for i in range(len(labels), len(bam_files)):
        labels.append('BAM%d' % (i+1))

    ############################################
    # intersect
    ############################################
    bam_fwd_coverages = []
    bam_rev_coverages = []
    for i in range(len(bam_files)):
        # make dir
        bam_dir = '%s/bam%d' % (options.out_dir,i+1)
        if not os.path.isdir(bam_dir):
            os.mkdir(bam_dir)

        # compute cov
        cov_fwd, cov_rev = compute_bam_cov(dfam_te, bam_files[i], bam_dir, options.stranded, options.min_cov, options.max_reads, options.msa_saved, options.clean)

        # smooth
        cov_fwd_smooth = kernel_smooth(cov_fwd)
        cov_rev_smooth = kernel_smooth(cov_rev)

        # add pseudocounts
        for i in range(len(cov_fwd_smooth)):
            cov_fwd_smooth[i] += 1
        for i in range(len(cov_rev_smooth)):
            cov_rev_smooth[i] += 1

        # save
        bam_fwd_coverages.append(cov_fwd_smooth)
        bam_rev_coverages.append(cov_rev_smooth)

    ############################################
    # plot
    ############################################
    print >> sys.stderr, 'Plotting.'
    
    plot_coverage(bam_fwd_coverages, labels, dfam_te, '%s/coverage_forward.pdf' % options.out_dir, reverse=False)
    plot_coverage(bam_rev_coverages, labels, dfam_te, '%s/coverage_reverse.pdf' % options.out_dir, reverse=True)

    # weblogo?
    '''
    if orient == 'fwd':
        subprocess.call('stockholm2fasta.py -c %s/hmmer_main/%s_%s.msa > %s/hmmer_main/%s_%s_cons.fa' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
    else:
        subprocess.call('stockholm2fasta.py -c %s/hmmer_main/%s_%s.msa > %s/hmmer_main/%s_%s_consr.fa' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
        subprocess.call('dna.py --rc %s/hmmer_main/%s_%s_consr.fa > %s/hmmer_main/%s_%s_cons.fa' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
        os.remove('%s/hmmer_main/%s_%s_consr.fa' % (options.out_dir,te,orient))

    subprocess.call('weblogo < %s/hmmer_main/%s_%s_cons.fa > %s/%s_%s_logo.eps' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
    os.remove('%s/hmmer_main/%s_%s_cons.fa' % (options.out_dir,te,orient))
    '''


################################################################################
# compute_bam_cov
#
# Compute coverage across TE seqeunces from a BAM file.
################################################################################
def compute_bam_cov(dfam_te, bam_file, bam_dir, stranded, min_cov, max_reads, msa_saved, clean):
    print >> sys.stderr, 'Computing %s coverage...' % bam_file

    if msa_saved:
        print >> sys.stderr, '\tNo %s renormalization with saved MSA.' % dfam_te
        renorm_fwd = 1.0
        renorm_rev = 1.0

    else:
        # hash reads by TE
        print >> sys.stderr, '\tHashing reads by %s.' % dfam_te
        te_reads = hash_te_reads(dfam_te, bam_file)

        # how long is the TE?
        hmm_file = '%s/hmms/%s.hmm' % (tempura.dfam_dir, dfam_te)
        for line in open(hmm_file):
            if line.startswith('LENG'):
                te_length = int(line[5:])
                break

        # do we see enough reads/len
        if len(te_reads)/float(te_length) < min_cov:
            renorm_fwd = 0
            renorm_rev = 0
        else:
            # make TE read fasta files
            print >> sys.stderr, '\tMaking %s read fasta files.' % dfam_te
            renorm_fwd, renorm_rev = make_te_read_fastas(bam_file, te_reads, bam_dir, stranded, max_reads)

    # compute coverage and re-normalize forward and reverse
    bam_cov_fwd = cov_renorm(dfam_te, '%s/fwd.fa' % bam_dir, renorm_fwd, msa_saved)
    if stranded:
        bam_cov_rev = cov_renorm(dfam_te, '%s/rev.fa' % bam_dir, renorm_rev, msa_saved)
    else:
        bam_cov_rev = []   # same as if we decided not to align it

    rm_if_is('%s/fwd.fa' % bam_dir)
    rm_if_is('%s/rev.fa' % bam_dir)
    if clean:
        rm_if_is('%s/fwd.msa' % bam_dir)
        rm_if_is('%s/rev.msa' % bam_dir)

    return bam_cov_fwd, bam_cov_rev


################################################################################
# cov_renorm
#
# Run hmmalign, extract coverage, and renormalize.
################################################################################
def cov_renorm(dfam_te, fasta_file, renorm, msa_saved):
    # hmmalign and parse
    bam_cov = hmmer_cov(dfam_te, fasta_file, msa_saved)

    # renormalize
    if len(bam_cov) > 0 and renorm != 1:
        print >> sys.stderr, '\tRenormalize %s coverage with %f' % (fasta_file, renorm)
        for i in range(len(bam_cov)):
            bam_cov[i] *= renorm

    return bam_cov


################################################################################
# hmmer_cov
#
# Run hmmalign on the fasta file versus the dfam_te profile HMM and compute
# coverage across the TE profile.
################################################################################
def hmmer_cov(dfam_te, fasta_file, msa_saved):
    bam_cov = []

    align_block = []
    block_cons_start = 0

    # hmmalign
    hmm_file = '%s/hmms/%s.hmm' % (tempura.dfam_dir, dfam_te)
    msa_file = '%s.msa' % os.path.splitext(fasta_file)[0]

    if not msa_saved and os.path.isfile(fasta_file):
        subprocess.call('hmmalign --dna %s %s > %s' % (hmm_file,fasta_file,msa_file), shell=True)
        
    if os.path.isfile(msa_file):
        msa_in = open(msa_file)

        line = msa_in.readline() # header
        for line in msa_in:
            if line.rstrip() in ['','//']:
                if align_block:
                    block_cons_start = process_align_block(align_block, bam_cov, block_cons_start)

                # reset align block
                align_block = []
            else:
                align_block.append(line)

        msa_in.close()

    return bam_cov


################################################################################
# kernel_smooth
#
# Ad hoc Gaussian kernel smoother.
################################################################################
def kernel_smooth(cov):
    # dnorm(-2:2, mean=0, sd=1)
    kernel_weights = [0.05, 0.25, 0.4, 0.25, 0.05]
    bw2 = (len(kernel_weights)-1)/2
    
    smooth_cov = []
    for i in range(len(cov)):
        smooth_cov.append(0)
        kernel_sum = 0
        for k in range(-bw2,bw2+1):
            if 0 <= i+k < len(cov):
                smooth_cov[-1] += cov[i+k]*kernel_weights[bw2+k]
                kernel_sum += kernel_weights[bw2+k]
        smooth_cov[-1] /= kernel_sum

    return smooth_cov


################################################################################
# hash_te_reads
#
# Use intersectBed to find the set of reads aligning to this TE in the BAM file,
# returning a dict mapping them to strand information.
#
# Input
#  dfam_te:
#  bam_file:
#  min_overlap:
#
# Output
#  te_reads: Dict mapping read header to a tuple of read strand and TE strand.
################################################################################
def hash_te_reads(dfam_te, bam_file, min_overlap=5):
    te_reads = {}

    dfam_gff = '%s/gffs/%s.gff' % (tempura.dfam_dir, dfam_te)

    p = subprocess.Popen('intersectBed -wo -bed -split -sorted -abam %s -b %s' % (bam_file, dfam_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')

        rheader = a[3]
        rstrand = a[5]

        tstrand = a[18]

        overlap = int(a[-1])

        # hash the overlap and the strands
        if overlap >= min_overlap:
            te_reads[rheader] = (rstrand,tstrand)
    p.communicate()

    return te_reads


################################################################################
# make_te_read_fastas
#
# Make a fasta file of overlapping reads for each TE.
#
# The strategy here is to loop through the BAM file, and print the reads
# as we encounter their alignments on the proper strand. Print only up to the
# max number of reads, but count everything, so we can't renormalize. Finally,
# filter the fasta files with too few reads.
################################################################################
def make_te_read_fastas(bam_file, te_reads, out_dir, stranded, max_reads):
    # open TE read fasta files
    fasta_fwd_out = open('%s/fwd.fa' % out_dir, 'w')
    fasta_rev_out = open('%s/rev.fa' % out_dir, 'w')

    # initialize counters for total reads
    total_fwd = 0
    total_rev = 0

    reads_printed = set()

    # print reads to fasta files
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        if aligned_read.qname in te_reads and (aligned_read.qname,aligned_read.is_read1) not in reads_printed:
            (rstrand, tstrand) = te_reads[aligned_read.qname]

            # only print if we match the read strand
            if (aligned_read.is_reverse and rstrand == '-') or (not aligned_read.is_reverse and rstrand == '+'):
                # TE determines reversal
                if tstrand == '+':
                    rseq = aligned_read.seq
                else:
                    rseq = dna.rc(aligned_read.seq)

                # count, and print
                if not stranded or rstrand == tstrand:
                    total_fwd += 1
                    if total_fwd < max_reads:
                        print >> fasta_fwd_out, '>%s\n%s' % (aligned_read.qname,rseq)
                else:
                    total_rev += 1
                    if total_rev < max_reads:
                        print >> fasta_rev_out, '>%s\n%s' % (aligned_read.qname,rseq)

                # save so we won't print again
                reads_printed.add((aligned_read.qname,aligned_read.is_read1))

    # close fasta files
    fasta_fwd_out.close()    
    fasta_rev_out.close()

    # return renormalization factors    
    renorm_fwd = max(1.0, total_fwd/float(max_reads))
    renorm_rev = max(1.0, total_rev/float(max_reads))

    return renorm_fwd, renorm_rev


################################################################################
# plot_coverage
#
# Plot normalized coverage across the TE for each of the BAM files.
################################################################################
def plot_coverage(bam_coverages, labels, dfam_te, out_pdf, reverse=False):
    df = {'indexes':[], 'coverage':[], 'data':[]}
    for i in range(len(bam_coverages)):
        if len(bam_coverages[i]) > 0:
            df['indexes'] += range(len(bam_coverages[i]))
            df['data'] += [labels[i]]*len(bam_coverages[i])

            cov_sum = float(sum(bam_coverages[i]))
            final_cov = [c/cov_sum for c in bam_coverages[i]]

            if reverse:
                df['coverage'] += final_cov[::-1]
            else:
                df['coverage'] += final_cov

    if len(df['indexes']) > 0:
        ggplot.plot('%s/te_cov.r' % tempura.r_dir, df, [dfam_te, out_pdf])


################################################################################
# process_align_block
#
# Process the block of multiple sequence alignment lines, incrementing the
# indexes in coverage and update the consensus start.
################################################################################
def process_align_block(align_block, coverage, block_cons_start):
    nts = ['A','C','G','T','a','c','g','t']

    # read consensus X markers
    consensus_x = align_block[-1].split()[2]

    # add zeroes to the coverage counters
    for cons_nt in consensus_x:
        if cons_nt == 'x':
            coverage.append(0)

    # for each block line
    for block_line in align_block:
        # if it's a true msa line
        if block_line[0] != '#':
            # get the msa
            block_a = block_line.split()
            rheader = block_a[0]
            rmsa = block_a[1]

            # rest the consensus index
            cons_i = block_cons_start

            # for each msa entry
            for msa_i in range(len(rmsa)):
                # if it's a consensus nt
                if consensus_x[msa_i] == 'x':
                    # if the msa has a valid entry
                    if rmsa[msa_i] in nts:
                        # increment coverage
                        coverage[cons_i] += 1
                    # increment the consensus index
                    cons_i += 1

    # update block consensus start to the pos reached
    return cons_i


################################################################################
# rm_if_is
#
# Check that the file exists before removing it.
################################################################################
def rm_if_is(f):
    if os.path.isfile(f):
        os.remove(f)


################################################################################
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
