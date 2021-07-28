#!/bin/env python3
####################################################################################################
## Copyright (C) 2016 Princess Margaret Bioinformatics and HPC Core - All Rights Reserved
## You may freely use, distribute and modify copies of this code within any of the systems currently 
## owned and operated by the University Health Network and the Bioinformatics and HPC Core. 
## If you require pieces of this code for any other purposes outside of these systems
## please contact the Bioinformatics and HPC Core directly for permission. 
##
## The Bioinformatics and HPC Core makes no claims as to the usability or reliability of 
## this code. Use for clinical purposes must only be done by those at UHN who are 
## authorized to do so using the Standard Operating Practices that have been developed specifically
## for this purpose.
####################################################################################################

"""
This script runs qualimap BAMQC on all the BAM files in a provided directory.
"""

import os
import argparse
import glob

import pmgctools
import qsub


def init():
    parser = argparse.ArgumentParser(
        description='This script runs qualimap BAMQC on all the BAM files in a provided directory.')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-f', '--format', default='pdf', help='Result format (pdf or html). Default: pdf')
    parser.add_argument('-g', '--bed', help='GTF/BED file of sequenced regions, needed for exome/targeted experiments')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-S', '--species', default="HUMAN", help='Species HUMAN or MOUSE default=HUMAN')
    parser.add_argument('-Q', '--queue', help='grid engine queue to use')
    args = parser.parse_args()
    return args


def bamqc(bam, outputdir, bed=None, species='HUMAN', format='PDF', waitlist=None, **other_qsub_options):
    """
    run bamqc on one bam file
    :param bam: bam file to QC
    :param outputdir: directory for the bamqc result
    :param bed: region file (GTF/BED)
    :param species: Species HUMAN or MOUSE
    :param mem: memory size
    :param format: Result format (pdf or html)
    :param other_qsub_options: other ooptions for qsub.qsub
    :return: job name
    """

    mem ='8G'

    modules = pmgctools.check_vars(['qualimap'])

    if not bam.endswith('.bam'):
        print('bamqc needs a bam file to process.')
        exit(1)
    sample = os.path.basename(bam)[:-4]

    # create output and tmp dir if not exists
    tmpdir = os.path.join(outputdir, 'tmp')
    sampledir = os.path.join(tmpdir, sample)
    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir, 0o770)

    if not bed:
        bed = pmgctools.get_var('bed')

    cmd = 'qualimap bamqc' if bed is None else 'qualimap bamqc -gff {}'.format(bed)
    cmd += ' -bam {} -gd {} -nt 4 -outdir {} -c --java-mem-size={}'.format(bam, species, sampledir, mem)
    if format.upper() == 'HTML':
        cmd += '\ncp -r {} {}\n'.format(sampledir, os.path.join(outputdir, sample + '_bamqc'))
    else:
        # need to run html version anyway for multi-sample bamqc
        cmd_pdf = cmd + ' -outformat PDF'
        cmd += '\n{}\nmv {} {}\n'.format(cmd_pdf, os.path.join(sampledir, 'report.pdf'),
                                         os.path.join(outputdir, sample + '_bamqc.pdf'))

    name = 'bamqc_' + sample
    wait = ",".join(waitlist) if waitlist is not None else None
    return qsub.qsub(name, cmd, modules=modules, waitlist=wait, other='cpu=4|mem={}|walltime=12:00:00'.format('16gb'), **other_qsub_options)


def multi_bamqc(samplelist, bamqc_waitlist, bamqcdir, format='PDF', **other_qsub_options):
    """
    run bamqc for multiple samples based on each bamqc result
    """
    modules = pmgctools.check_vars(['qualimap'])

    # create sample list
    names = []
    samplelist_file = os.path.join(bamqcdir, 'samplelist')
    with open(samplelist_file, 'w') as f:
        for sample in samplelist:
            names.append(sample)
            print('{}\t{}'.format(sample, os.path.join(bamqcdir, 'tmp', sample)), file=f)

    name = '__'.join(names)[:200]
    multi_dir = os.path.join(bamqcdir, 'tmp', 'multi_bamqc')
    cmd = "qualimap multi-bamqc -d {} -outdir {}".format(samplelist_file, multi_dir)
    if format.upper() == 'HTML':
        cmd += '\nmv {} {}\n'.format(multi_dir, os.path.join(bamqcdir, name + '_bamqc'))
    else:
        cmd += ' -outformat pdf'
        cmd += '\nmv {} {}\n'.format(os.path.join(multi_dir, 'report.pdf'),
                                     os.path.join(bamqcdir, name + '_bamqc.pdf'))

    cmd += '''while read line; do sample=( $line ); echo ${{sample[0]}}>>{output}/coverage.txt; echo ---->>{output}/coverage.txt; grep 'mean coverageData\|std coverageData' ${{sample[1]}}/genome_results.txt>>{output}/coverage.txt; done<{output}/samplelist; rm -rf {output}/tmp; rm {output}/samplelist'''.format(output=bamqcdir)

    return qsub.qsub('bamqc_' + name, cmd, modules=modules, waitlist=','.join(bamqc_waitlist), other='cpu=1|mem={}|walltime=2:00:00'.format('16gb'), **other_qsub_options)


if __name__ == '__main__':
    args = init()
    if args.ini:
        pmgctools.read_vars(args.ini)


    waitlist = []
    samplelist = []
    for sample in glob.glob(os.path.join(args.source, '*.bam')):
        samplelist.append(os.path.basename(sample)[:-4])
        waitlist.append(bamqc(sample, args.output, bed=args.bed, species=args.species, format=args.format,
                              log=args.log, qsub=args.qsub, dry=args.dry, queue=args.queue))

    multi_bamqc(samplelist, waitlist, bamqcdir=args.output, format=args.format, log=args.log, qsub=args.qsub, dry=args.dry, queue=args.queue)
