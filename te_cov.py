#!/usr/bin/env python
from optparse import OptionParser
import copy, glob, os, random, re, shutil, stats, subprocess, sys, pdb, tempfile
import pysam
import bam_fragments, dna, gff, ggplot, te, util
import tempura

################################################################################
# te_cov.py
#
# Plot BAM coverage across a TE family.
#
# TODO:
#  (1) Fix paired end reads w/ same qname.
#      (a) It's a mess. You can't hash which read it is in the intersectBed
#           because you have to use -split -bed. But if the paired reads have
#           the same qname, it will get confused about the orientation they
#           should be in.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <bam1> ... <bam2> ... <bamN>'
    parser = OptionParser(usage)

    parser.add_option('-l', dest='labels', help='BAM labels (comma separated) for plots')
    parser.add_option('-s', dest='stranded', default=False, action='store_true', help='Seq is stranded so consider orientation [Default: %default]')
    parser.add_option('-t', dest='te_gff', default='%s/te.gff'%tempura.dfam_dir, help='TE GFF (use a single TE GFF if desired) [Default: %default]')

    # MSA
    parser.add_option('-c', dest='clean', default=False, action='store_true', help='Clean up removing temp files like the MSA computed for the reads [Default: %default]')
    parser.add_option('-m', dest='msa_saved', default=False, action='store_true', help='MSAs computed by HMMer are already done and in the output directory [Default: %default]')

    # limitations
    parser.add_option('-n', dest='num_plots', default=100, type='int', help='Number of plots to make [Default: %default]')
    parser.add_option('-a', dest='max_reads', default=10000, type='int', help='Maximum number of reads to hmmalign [Default: %default]')

    parser.add_option('-o', dest='out_dir', default='te_cov', help='Output directory [Default: %default]')

    (options,args) = parser.parse_args()

    if len(args) < 2:
        parser.error('Must provide BAM file(s)')
    else:
        bam_files = args

    ############################################
    # prep
    ############################################
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    if options.labels:
        labels = options.labels.split(',')
    else:
        labels = []

    for i in range(len(labels), len(bam_files)):
        labels.append('BAM%d' % (i+1))

    ############################################
    # intersect
    ############################################
    bam_te_coverages = []
    for i in range(len(bam_files)):
        # make dir
        bam_dir = '%s/bam%d' % (options.out_dir,i+1)
        if not os.path.isdir(bam_dir):
            os.mkdir(bam_dir)

        # compute coverage
        te_bam_cov = compute_bam_cov(options.te_gff, bam_files[i], bam_dir, options.stranded, options.max_reads, options.msa_saved, options.clean)

        # process coverage
        bam_frags = float(bam_fragments.count(bam_files[i]))
        for te_key in te_bam_cov:
            # smooth
            te_bam_cov[te_key] = kernel_smooth(te_bam_cov[te_key])

            # add pseudocounts (why?)
            #for i in range(len(te_cov_smooth[te_key])):
            #    te_cov_smooth[te_key][i] += 1

            # normalize            
            for c in range(len(te_bam_cov[te_key])):
                te_bam_cov[te_key][c] /= bam_frags

        # save
        bam_te_coverages.append(te_bam_cov)

    ############################################
    # choose top TEs to plot
    ############################################
    te_plot_scores = []
    
    all_tes = set()
    for i in range(len(bam_te_coverages)):
        for te_key in bam_te_coverages[i]:
            all_tes.add(te_key)

    # this must have nothing, right?
    print 'len(all_tes) = %d' % len(all_tes)

    for te_key in all_tes:
        plot_score = max([0]+[max([0]+bam_te_coverages[i].get(te_key,[])) for i in range(len(bam_te_coverages))])
        te_plot_scores.append((plot_score, te_key))

    te_plot_scores.sort(reverse=True)

    plot_tes = [te_key for (plot_score, te_key) in te_plot_scores[:options.num_plots]]


    ############################################
    # plot
    ############################################
    print >> sys.stderr, 'Plotting.'
    
    for dfam_te, orient in plot_tes:
        print >> sys.stderr, '\t%s %s' % (dfam_te, orient)
        plot_coverage(bam_te_coverages, dfam_te, orient, labels, options.out_dir)

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
def compute_bam_cov(te_gff, bam_file, bam_dir, stranded, max_reads, msa_saved, clean):
    print >> sys.stderr, 'Computing %s coverage...' % bam_file

    if msa_saved:
        print >> sys.stderr, '\tNo renormalization with saved MSA.'

        te_renorm = {}

        te_re = re.compile('%s/(.+)_(fwd|rev)\.msa' % bam_dir)
        for msa_file in glob.glob('%s/*_*.msa' % bam_dir):
            te_match = te_re.search(msa_file)

            te = te_match.group(1)
            orient = te_match.group(2)

            te_renorm[(te,orient)] = 1.0

    else:
        # hash reads by TE
        print >> sys.stderr, '\tHashing reads by repeat.'
        read_tes = hash_te_reads(te_gff, bam_file)

        # make TE read fasta files
        print >> sys.stderr, '\tMaking read fasta files.'
        te_renorm = make_te_read_fastas(te_gff, bam_file, read_tes, bam_dir, stranded, max_reads)

    te_bam_cov = {}
    for dfam_te, orient in te_renorm:
        fasta_file = '%s/%s_%s.fa' % (bam_dir, dfam_te, orient)

        # compute coverage and re-normalize
        te_bam_cov[(dfam_te,orient)] = cov_renorm(dfam_te, fasta_file, te_renorm[(dfam_te,orient)], msa_saved)

        # consider cleaning
        rm_if_is(fasta_file)
        if clean:
            rm_if_is('%s/%s_%s.msa' % (bam_dir, dfam_te, orient))

    return te_bam_cov


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
#  bam_file:
#  te_gff:
#  min_overlap:
#
# Output
#  te_reads:     Dict mapping read header to a dict mapping repeats to tuples of
#                 read strand and TE strand.
################################################################################
def hash_te_reads(te_gff, bam_file, min_overlap=5):
    te_reads = {}
    
    p = subprocess.Popen('intersectBed -wo -bed -split -sorted -abam %s -b %s' % (bam_file, te_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')

        rheader = a[3]
        rstrand = a[5]

        tstrand = a[18]
        dfam_te = gff.gtf_kv(a[20])['dfam']

        overlap = int(a[-1])

        # hash the overlap and the strands
        if overlap >= min_overlap:
            te_reads.setdefault(rheader,{})[dfam_te] = (rstrand, tstrand)
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
#
# If stranded==False, just use the 'fwd' file.
#
# Output
#  te_renorm: dict mapping (dfam_te,orient) tuples to a renormalization factor.
################################################################################
def make_te_read_fastas(te_gff, bam_file, read_tes, out_dir, stranded, max_reads):
    # open TE read fasta files
    te_fastas = {}
    for line in open(te_gff):
        a = line.split('\t')
        dfam_te = gff.gtf_kv(a[8])['dfam']
        if not (dfam_te,'fwd') in te_fastas:
            te_fastas[(dfam_te,'fwd')] = open('%s/%s_fwd.fa' % (out_dir,dfam_te), 'w')
            te_fastas[(dfam_te,'rev')] = open('%s/%s_rev.fa' % (out_dir,dfam_te), 'w')

    # initialize counters for total reads
    te_totals = {}
    for dfam_te, orient in te_fastas:
        te_totals[dfam_te, orient] = 0

    # print reads to fasta files
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        this_read_tes = read_tes.get(aligned_read.qname,{})

        for dfam_te in this_read_tes.keys():
            if this_read_tes[dfam_te] != None:
                (rstrand, tstrand) = this_read_tes[dfam_te]

                # only print if we match the read strand
                if (aligned_read.is_reverse and rstrand == '-') or (not aligned_read.is_reverse and rstrand == '+'):
                    # TE determines reversal
                    if tstrand == '+':
                        rseq = aligned_read.seq
                    else:
                        rseq = dna.rc(aligned_read.seq)

                    # count, and print
                    if not stranded or rstrand == tstrand:
                        te_totals[(dfam_te,'fwd')] += 1
                        if te_totals[(dfam_te,'fwd')] < max_reads:
                            print >> te_fastas[(dfam_te,'fwd')], '>%s\n%s' % (aligned_read.qname,rseq)
                    else:
                        te_totals[(dfam_te,'rev')] += 1
                        if te_totals[(dfam_te,'rev')] < max_reads:
                            print >> te_fastas[(dfam_te,'rev')], '>%s\n%s' % (aligned_read.qname,rseq)

                    # specify printed
                    this_read_tes[dfam_te] = None

    # post-process fasta files
    te_renorm = {}
    for dfam_te, orient in te_fastas:
        # close
        te_fastas[(dfam_te, orient)].close()

        # return renormalization factors
        if te_totals[(dfam_te,orient)] > 10:
            te_renorm[(dfam_te,orient)] = max(1.0, te_totals[(dfam_te,orient)]/float(max_reads))

    return te_renorm


################################################################################
# plot_coverage
#
# Plot normalized coverage across the TE for each of the BAM files.
################################################################################
def plot_coverage(bam_te_coverages, dfam_te, orient, labels, out_dir):

    df = {'indexes':[], 'coverage':[], 'coverage_norm':[], 'data':[]}
    for i in range(len(bam_te_coverages)):        
        bam_coverage = bam_te_coverages[i].get((dfam_te,orient),[])

        if len(bam_coverage) > 0:
            df['indexes'] += range(len(bam_coverage))
            df['data'] += [labels[i]]*len(bam_coverage)

            cov_sum = float(sum(bam_coverage))
            bam_coverage_norm = [c/cov_sum for c in bam_coverage]

            if orient == 'rev':
                df['coverage'] += bam_coverage[::-1]
                df['coverage_norm'] += bam_coverage_norm[::-1]
            else:
                df['coverage'] += bam_coverage
                df['coverage_norm'] += bam_coverage_norm

    if len(df['indexes']) > 0:
        out_pre = '%s/%s_%s_cov' % (out_dir, dfam_te, orient)
        ggplot.plot('%s/te_cov.r' % tempura.r_dir, df, [dfam_te, out_pre])


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
