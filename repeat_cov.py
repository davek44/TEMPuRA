#!/usr/bin/env python
from optparse import OptionParser
import glob, os, random, re, shutil, stats, subprocess, sys, pdb, tempfile
import pysam
import dna, gff, ggplot, te

################################################################################
# te_bam_cov_sim.py
#
# Plot BAM coverage across TE families, normalized by simulated reads.
#
# There's an annoying bug here where multiple RepeatMasker repeats can be
# mapped to the same DFAM repeat (e.g. THE1A and THE1A-int to THE1A) and the
# program gets confused downstream. The solution (that I didn't want to implement
# this late in the game) is to hash reads to DFAM TEs instead of RM TEs. Then
# it'll be easy to print to DFAM-named FASTA files and MSA files, and things
# flow from there.
################################################################################

################################################################################
# main
################################################################################
def main():
    usage = 'usage: %prog [options] <feature_bam> <control_bam>'
    parser = OptionParser(usage)
    parser.add_option('-m', dest='msa_done', default=False, action='store_true', help='MSAs computed by HMMer are already done and in the output directory [Default: %default]')
    parser.add_option('-o', dest='out_dir', default='cov', help='Output directory [Default: %default]')
    parser.add_option('-p', dest='plot_num', default=300, type='int', help='Number of plots to make [Default: %default]')
    parser.add_option('-r', dest='te_gff', default='%s/hg19.fa.out.tp.sort.gff' % os.environ['MASK'])
    (options,args) = parser.parse_args()

    if len(args) != 2:
        parser.error('Must provide feature BAM and control BAM')
    else:
        feature_bam = args[0]
        control_bam = args[1]

    ############################################
    # prepare output dir
    ############################################
    if not os.path.isdir(options.out_dir):
        os.mkdir(options.out_dir)

    ############################################
    # intersect
    ############################################
    te_feature_cov = compute_bam_cov(feature_bam, options.te_gff, '%s/hmmer_main'%options.out_dir, options.msa_done)
    te_null_cov = compute_bam_cov(control_bam, options.te_gff, '%s/hmmer_null'%options.out_dir, options.msa_done)

    # add null pseudocounts
    for te_ornt in te_feature_cov:
        te_len = len(te_feature_cov[te_ornt])

        if not te_ornt in te_null_cov:
            te_null_cov[te_ornt] = [0]*te_len

        for i in range(te_len):
            te_null_cov[te_ornt][i] += 1

    # smooth
    for te_ornt in te_feature_cov:
        te_feature_cov[te_ornt] = kernel_smooth(te_feature_cov[te_ornt])
    for te_ornt in te_null_cov:
        te_null_cov[te_ornt] = kernel_smooth(te_null_cov[te_ornt])

    ############################################
    # plot
    ############################################
    # sort by max coverage and do the top ones
    #te_mean_coverages = [(stats.mean(te_feature_cov[te_ornt]),te_ornt) for te_ornt in te_feature_cov]
    #te_mean_coverages.sort(reverse=True)
    #plot_tes = [te_ornt for (te_mean_cov,te_ornt) in te_mean_coverages[:options.plot_num]]

    '''
    for te_ornt in te_feature_cov:
        te, orient = te_ornt
        print >> sys.stderr, te, orient, te_feature_cov[te_ornt]
    '''

    te_max_coverages = [(max([1.0]+te_feature_cov[te_ornt]),te_ornt) for te_ornt in te_feature_cov]
    te_max_coverages.sort(reverse=True)
    plot_tes = [te_ornt for (te_max_cov,te_ornt) in te_max_coverages[:options.plot_num] if te_max_cov > 200]

    print >> sys.stderr, 'Determined %d TEs to plot' % len(plot_tes)

    for te_ornt in plot_tes:
        te, orient = te_ornt

        print >> sys.stderr, 'Plotting %s %s...' % (te,orient),

        te_len = len(te_feature_cov[te_ornt])

        if orient == 'fwd':
            tfc = te_feature_cov[te_ornt]
            tnc = te_null_cov[te_ornt]
        else:
            tfc = te_feature_cov[te_ornt][::-1]
            tnc = te_null_cov[te_ornt][::-1]

        # BAM coverage
        df = {'indexes':range(te_len), 'coverage':tfc}
        out_pdf = '%s/%s_%s_main.pdf' % (options.out_dir,te.replace('/','-'),orient)
        ggplot.plot('%s/te_bam_coverage_bam.r'%os.environ['RDIR'], df, [te,out_pdf])

        # sim coverage
        df = {'indexes':range(te_len), 'coverage':tnc}
        out_pdf = '%s/%s_%s_null.pdf' % (options.out_dir,te.replace('/','-'),orient)
        ggplot.plot('%s/te_bam_coverage_bam.r'%os.environ['RDIR'], df, [te,out_pdf])

        # BAM normalized coverage
        feature_sum = float(sum(tfc))
        feature_norm = [tfc[i]/feature_sum for i in range(te_len)]
        null_sum = float(sum(tnc))
        null_norm = [tnc[i]/null_sum for i in range(te_len)]

        norm_cov = [feature_norm[i] - null_norm[i] for i in range(te_len)]
        df = {'indexes':range(te_len), 'coverage':norm_cov}
        out_pdf = '%s/%s_%s_norm.pdf' % (options.out_dir,te.replace('/','-'),orient)
        ggplot.plot('%s/te_bam_coverage_bam.r'%os.environ['RDIR'], df, [te,out_pdf])

        # weblogo
        if orient == 'fwd':
            subprocess.call('stockholm2fasta.py -c %s/hmmer_main/%s_%s.msa > %s/hmmer_main/%s_%s_cons.fa' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
        else:
            subprocess.call('stockholm2fasta.py -c %s/hmmer_main/%s_%s.msa > %s/hmmer_main/%s_%s_consr.fa' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
            subprocess.call('dna.py --rc %s/hmmer_main/%s_%s_consr.fa > %s/hmmer_main/%s_%s_cons.fa' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
            os.remove('%s/hmmer_main/%s_%s_consr.fa' % (options.out_dir,te,orient))

        subprocess.call('weblogo < %s/hmmer_main/%s_%s_cons.fa > %s/%s_%s_logo.eps' % (options.out_dir,te,orient,options.out_dir,te,orient), shell=True)
        os.remove('%s/hmmer_main/%s_%s_cons.fa' % (options.out_dir,te,orient))

        print >> sys.stderr, 'Done'


################################################################################
# compute_bam_cov
#
# Compute coverage across TE seqeunces from a BAM file.
################################################################################
def compute_bam_cov(bam_file, te_gff, hmmer_dir=None, msa_done=False):
    print >> sys.stderr, 'Computing %s coverage...' % bam_file

    hmmer_temp = False

    if not msa_done:
        # create temp directory        
        if hmmer_dir == None:
            hmmer_dir = tempfile.mkdtemp()
            hmmer_temp = True
        else:
            if not os.path.isdir(hmmer_dir):
                os.mkdir(hmmer_dir)

        # hash reads by TE
        print >> sys.stderr, '\tHashing reads by TE...',
        read_tes = hash_te_reads(bam_file, te_gff)
        print >> sys.stderr, 'Done'

        # make TE read fasta files
        print >> sys.stderr, '\tMaking TE read fasta files...',
        te_renorm = make_te_read_fastas(bam_file, te_gff, read_tes, hmmer_dir, read_max=30000)
        print >> sys.stderr, 'Done'

    else:        
        # find already done combinations and renormalization factor
        te_renorm = msa_renorm(bam_file, te_gff, hmmer_dir)

        '''
        # missing renorm sucks, but this is also wrong
        # because it stores dfam names, whereas above
        # stores repeatmasker

        print >> sys.stderr, 'Warning: TE fragment count renormalizations are missing'

        te_re = re.compile('%s/(.+)_(fwd|rev)\.msa' % hmmer_dir)
        te_renorm = {}
        for msa_file in glob.glob('%s/*_*.msa' % hmmer_dir):
            te_match = te_re.search(msa_file)

            te_repeat = te_match.group(1)
            orient = te_match.group(2)

            te_renorm[(te_repeat,orient)] = 1.0
        '''

    # hmmalign and parse
    te_bam_cov = {}    
    for te_repeat, orient in te_renorm:
        print >> sys.stderr, 'hmmalign %s %s reads...' % (te_repeat, orient),

        fasta_file = '%s/%s_%s.fa' % (hmmer_dir,te_repeat,orient)

        # map RepeatMasker to DFAM
        dfam_tes = map_rm_dfam(te_repeat)
        for dfam_te in dfam_tes:
            
            # because of renorm and weird names, gotta check
            msa_file = '%s/%s_%s.msa' % (hmmer_dir,dfam_te,orient)
            if os.path.isfile(msa_file):

                # run HMMer and compute coverage
                hmmer_cov(fasta_file, dfam_te, orient, te_bam_cov, hmmer_dir, msa_done)

                if te_renorm[(te_repeat,orient)] != 1.0:
                    # renormalize
                    print >> sys.stderr, ' renormalize with %f' % te_renorm[(te_repeat,orient)]
                    for i in range(len(te_bam_cov[(dfam_te,orient)])):
                        te_bam_cov[(dfam_te,orient)][i] *= te_renorm[(te_repeat,orient)]

                # clean as we go
                if hmmer_temp:
                    msa_file = '%s/%s_%s.msa' % (hmmer_dir,dfam_te,orient)
                    os.remove(msa_file)

        # clean as we go
        if not msa_done:
            os.remove(fasta_file)

        print >> sys.stderr, 'Done'

    # clean
    if hmmer_temp:
        shutil.rmtree(hmmer_dir)

    return te_bam_cov


################################################################################
# count_msa_reads
#
# Count the number of reads in a MSA file.
################################################################################
def count_msa_reads(msa_file):
    reads = 0

    msa_in = open(msa_file)
    line = msa_in.readline() # header
    line = msa_in.readline() # blank
    line = msa_in.readline()

    while line and line.rstrip() not in ['','//']:
        if line[0] != '#':
            reads += 1
        line = msa_in.readline()

    msa_in.close()

    return reads


################################################################################
# estimate_read_stats
#
# Compute mean read length by sampling the first N reads.
################################################################################
def estimate_read_stats(bam_file):
    samples = 100000
    s = 0
    read_lengths = []
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        read_lengths.append(aligned_read.rlen)
        s += 1
        if s >= samples:
            break
            
    mean_f, sd_f = stats.mean_sd(read_lengths)

    return int(mean_f+0.5), int(sd_f+0.5)


################################################################################
# hmmer_cov
#
# Run hmmalign on the fasta file versus the dfam_te profile HMM and compute
# coverage across the TE profile.
################################################################################
def hmmer_cov(fasta_file, dfam_te, orient, te_bam_cov, hmmer_dir, msa_done):
    te_bam_cov[(dfam_te,orient)] = []

    align_block = []
    block_cons_start = 0

    # hmmalign
    hmm_file = '%s/hmms/%s.hmm' % (os.environ['DFAM'],dfam_te)
    msa_file = '%s/%s_%s.msa' % (hmmer_dir,dfam_te,orient)

    if not msa_done:
        subprocess.call('hmmalign --dna %s %s > %s' % (hmm_file,fasta_file,msa_file), shell=True)
        
    msa_in = open(msa_file)
    line = msa_in.readline() # header
    for line in msa_in:
        if line.rstrip() in ['','//']:
            if align_block:
                block_cons_start = process_align_block(align_block, te_bam_cov[(dfam_te,orient)], block_cons_start)

            # reset align block
            align_block = []
        else:
            align_block.append(line)

    msa_in.close()


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
# map_rm_dfam
#
# Map a RepeatMasker name to a DFAM name.
################################################################################
def map_rm_dfam(repeat):
    if os.path.isfile('%s/hmms/%s.hmm' % (os.environ['DFAM'],repeat)):
        dfam_reps = [repeat]
    elif os.path.isfile('%s/hmms/%sv.hmm' % (os.environ['DFAM'],repeat)):
        dfam_reps = [repeat+'v']
    else:
        hmm_files = glob.glob('%s/hmms/%s_*.hmm' % (os.environ['DFAM'],repeat))

        # if no hits
        if len(hmm_files) == 0:
            # try removing "-int"
            if repeat[-4:] == '-int' and os.path.isfile('%s/hmms/%s.hmm' % (os.environ['DFAM'],repeat[:-4])):
                dfam_reps = [repeat[:-4]]
            else:
                # missing
                print >> sys.stderr, 'Missing DFAM name for %s' % repeat
                dfam_reps = []

        # if hits
        else:
            # grab em
            dfam_reps = []
            for i in range(len(hmm_files)):
                start = hmm_files[i].rfind('/')+1
                end = hmm_files[i].rfind('.hmm')
                dfam_reps.append(hmm_files[i][start:end])

    return dfam_reps


################################################################################
# hash_te_reads
#
# Return a dict keyed by read headers pointing to another dict that maps
# TE repeat names to the respective strands of the overlap
################################################################################
def hash_te_reads(bam_file, te_gff, min_overlap=5):
    read_tes = {}

    #p = subprocess.Popen('intersectBed -wo -bed -split -abam %s -b %s' % (bam_file, te_gff), shell=True, stdout=subprocess.PIPE)
    p = subprocess.Popen('intersectBed -wo -bed -split -sorted -abam %s -b %s' % (bam_file, te_gff), shell=True, stdout=subprocess.PIPE)
    for line in p.stdout:
        a = line.split('\t')

        rheader = a[3]
        rstrand = a[5]

        tstrand = a[18]
        te_repeat = gff.gtf_kv(a[20])['repeat']

        overlap = int(a[-1])

        # hash the overlap and the strands
        if overlap >= min_overlap:
            read_tes.setdefault(rheader,{})[te_repeat] = (rstrand,tstrand)
    p.communicate()

    return read_tes


################################################################################
# make_te_read_fastas
#
# Make a fasta file of overlapping reads for each TE.
################################################################################
def make_te_read_fastas(bam_file, te_gff, read_tes, out_dir, read_max=30000, read_min=100):
    # count reads per TE
    te_reads = {}
    for read_id in read_tes:
        for te in read_tes[read_id]:
            te_reads[te] = te_reads.get(te,0) + 1

    # compute a sampling probability for each TE
    te_sample = {}
    for te in te_reads:
        te_sample[te] = min(1.0, read_max / float(te_reads[te]))

    # open TE read fasta files
    te_fastas = {}
    for line in open(te_gff):
        a = line.split('\t')
        te_repeat = gff.gtf_kv(a[8])['repeat']
        if not (te_repeat,'fwd') in te_fastas:
            te_fastas[(te_repeat,'fwd')] = open('%s/%s_fwd.fa' % (out_dir,te_repeat), 'w')
            te_fastas[(te_repeat,'rev')] = open('%s/%s_rev.fa' % (out_dir,te_repeat), 'w')

    # make read fasta files
    te_reads_all = {}
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        this_read_tes = read_tes.get(aligned_read.qname,{})
        for te_repeat in this_read_tes.keys():
            # if not yet printed
            if this_read_tes[te_repeat] != None:
                (rstrand, tstrand) = this_read_tes[te_repeat]

                # only print if we match the read strand
                if (aligned_read.is_reverse and rstrand == '-') or (not aligned_read.is_reverse and rstrand == '+'):
                    # TE determines reversal
                    if tstrand == '+':
                        rseq = aligned_read.seq
                    else:
                        rseq = dna.rc(aligned_read.seq)

                    # count
                    if rstrand == tstrand:
                        te_reads_all[(te_repeat,'fwd')] = te_reads_all.get((te_repeat,'fwd'),0) + 1
                    else:
                        te_reads_all[(te_repeat,'rev')] = te_reads_all.get((te_repeat,'rev'),0) + 1

                    # maybe print
                    if random.random() < te_sample[te_repeat]:
                        if rstrand == tstrand:
                            print >> te_fastas[(te_repeat,'fwd')], '>%s\n%s' % (aligned_read.qname,rseq)
                        else:
                            print >> te_fastas[(te_repeat,'rev')], '>%s\n%s' % (aligned_read.qname,rseq)

                    # specify printed
                    this_read_tes[te_repeat] = None
    
    # close fasta files, trim near-empties, and hash renormalization factors
    te_renorm = {}
    for te_ornt in te_fastas:
        # close
        te_fastas[te_ornt].close()

        fasta_file = '%s/%s_%s.fa' % (out_dir,te_ornt[0],te_ornt[1])

        # count reads
        num_reads = 0
        for line in open(fasta_file):
            if line[0] == '>':
                num_reads += 1

        # remove if too few
        if num_reads < read_min:
            os.remove(fasta_file)
        else:
            te_renorm[te_ornt] = te_reads_all[te_ornt] / float(num_reads)

    return te_renorm


################################################################################
# msa_renorm
#
# Compute the TE renormalization factors for a finished MSA.
################################################################################
def msa_renorm(bam_file, te_gff, hmmer_dir):
    # map reads to TEs
    read_tes = hash_te_reads(bam_file, te_gff)
    
    # count reads in BAM
    te_reads_all = {}
    for aligned_read in pysam.Samfile(bam_file, 'rb'):
        this_read_tes = read_tes.get(aligned_read.qname,{})
        for te_repeat in this_read_tes.keys():
            # if not yet counted
            if this_read_tes[te_repeat] != None:
                (rstrand, tstrand) = this_read_tes[te_repeat]

                # only count if we match the read strand
                if (aligned_read.is_reverse and rstrand == '-') or (not aligned_read.is_reverse and rstrand == '+'):
                    # count
                    if rstrand == tstrand:
                        te_reads_all[(te_repeat,'fwd')] = te_reads_all.get((te_repeat,'fwd'),0) + 1
                    else:
                        te_reads_all[(te_repeat,'rev')] = te_reads_all.get((te_repeat,'rev'),0) + 1

                    # specify counted
                    this_read_tes[te_repeat] = None

    # count reads in MSAs
    te_re = re.compile('%s/(.+)_(fwd|rev)\.msa' % hmmer_dir)
    te_reads_msa = {}
    for msa_file in glob.glob('%s/*_*.msa' % hmmer_dir):
        te_match = te_re.search(msa_file)

        dfam_te = te_match.group(1)
        orient = te_match.group(2)

        te_reads_msa[(dfam_te,orient)] = count_msa_reads(msa_file)


    ##################################################
    # Ugh, needed this to avoid nasty problems with
    # two repeats that map to the same DFAM
    ##################################################
    buggy_repeats = set()
    dfam_repeat = te.map_dfam_repeat()
    for dfam in dfam_repeat:
        if len(dfam_repeat[dfam]) > 1:
            for repeat in dfam_repeat[dfam]:
                buggy_repeats.add(repeat)
            

    # compute renormalization factors
    te_renorm = {}
    for te_ornt in te_reads_all:
        te_repeat, orient = te_ornt

        if te_repeat in buggy_repeats:
            te_renorm[te_ornt] = 1.0
            print >> sys.stderr, te_ornt, 'buggy'
        
        else:
            # map to dfam id and check for consistency
            reads_msa = None
            dfam_tes = map_rm_dfam(te_repeat)
            for dfam_te in dfam_tes:
                reads_msa_dfam = te_reads_msa.get((dfam_te, orient),0)

                if reads_msa == None:
                    reads_msa = reads_msa_dfam
                else:
                    if reads_msa != reads_msa_dfam:
                        print >> sys.stderr, '%s/%s differing msa counts %d vs %d' % (te_repeat, orient, reads_msa, reads_msa_dfam)
                        exit(1)

            # if found, save
            if reads_msa != None and reads_msa > 0:
                te_renorm[te_ornt] = te_reads_all[te_ornt] / float(reads_msa)
                print >> sys.stderr, te_ornt, te_reads_all[te_ornt], reads_msa

    return te_renorm


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
# __main__
################################################################################
if __name__ == '__main__':
    main()
    #pdb.runcall(main)
