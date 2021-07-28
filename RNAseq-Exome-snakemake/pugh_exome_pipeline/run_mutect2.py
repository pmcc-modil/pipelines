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
This script runs mutect on tumor/normal pair of bam files
"""

import sys
import os
import argparse
import re

import qsub
import pmgctools


def init():
    parser = argparse.ArgumentParser(
        description='This script runs mutect on tumor/normal pair of bam files')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit to.')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-d', '--dbsnp',
                        help='current dbSNP reference in VCF format (if sample on multiple lanes and need to merge)')
    parser.add_argument('-m', '--cosmic',
                        help=' COSMIC vcf file)')
    parser.add_argument('-c', '--coverage',
                        help="coverage to downsample to at a given locus, enter an integer (e.g. 1000) or gatk to "
                             "use GATK default settings. If not provided, no downsampling will be performed.")
    parser.add_argument('-g', '--bed',
                        help='BED file containing intervals of interest'
                        )
    parser.add_argument('--no-padding', dest='no_padding', action='store_true', default=False,
                        help='We use default padding of 100 for intervals, this option will disable padding.'
                        )
    parser.add_argument('-C', '--config', required=True,
                        help='mutect config file. Put tumor and normal sample pair in one line seperated by space or , '
                             'i.e tumor1.bam,normal1.bam')
    parser.add_argument('-a', '--args', default='', help='additional arguments to be provided to mutect')

    options = parser.parse_args()
    return options


def mutect(source, outputdir, sample, bed=None, no_padding=False, dbsnp=None, cosmic=None, coverage=None, args='', waitlist=None, extension=None,
           **other_qsub_options):
    """
    run mutect on tumor, normal pair samples
    :param source: directory contains bam files
    :param outputdir: output directory
    :param sample: comma (,) separated tumor sample pair or tumor sample only. tumor sample first
    :param bed: region file
    :param dbsnp: dbsno vcf file
    :param cosmic: comsic vcf file
    :param coverage: downsampling see source code
    :param args: additional arguments to be provided to mutect
    :param other_qsub_options: other options for qsub script
    :return: job name, mutect_<sample>
    """
    tools = ['gatk', 'vcftools']
    modules = pmgctools.check_vars(tools)
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())
    wait = ",".join(waitlist) if waitlist is not None else None

    pair = sample.split(",")
    if not extension:
        if os.path.exists(os.path.join(source, pair[0] + '.processed.bam')):
            extension = '.processed.bam'
        elif os.path.exists(os.path.join(source, pair[0] + '.bam')):
            extension = '.bam'
        else:
            print('ERROR: I can not find bam files for samples "{}"'.format(sample), file=sys.stderr)
            exit(1)

    tumor_bam = os.path.join(source, pair[0] + extension)
    input = "-I:tumor {} ".format(tumor_bam)
    if len(pair) == 2:
        normal_bam = os.path.join(source, pair[1] + extension)
        input += "-I:normal {} ".format(normal_bam)
    elif len(pair) == 1:
        print("WARNING: only one sample found. I am going to use it as tumor sample.", file=sys.stderr)
        normal_bam = None
    else:
        print("ERROR: Multiple samples found. Please check the conf file.", file=sys.stderr)
        exit(1)


    if not bed:
        bed = pmgctools.get_var('bed')
    if bed:
        args += " --intervals " + bed + " --interval_padding "
        args += '0' if no_padding else '100'

    if not coverage:
        coverage = "-dt None"
    elif coverage == "gatk":
        coverage = ""
    elif coverage.isdigit():
        coverage = "-dcov " + coverage
    else:
        print("Error: wrong downsampling value " + coverage)

    ref = pmgctools.check_var("REF")

    if not cosmic:
        cosmic = pmgctools.get_var("COSMIC")
    cosmic = '--cosmic ' + cosmic if cosmic else ''

    if not dbsnp:
        dbsnp = pmgctools.check_var("dbSNP")

    sample = sample.replace(',', '__')[:200]
#    output = os.path.join(outputdir, sample + '.mutect.out')
    vcf = os.path.join(outputdir, sample + '.mutect.vcf')
#    coverage_file = os.path.join(outputdir, sample + '.coverage.wig.txt')
#    filtered = os.path.join(outputdir, sample + '.mutect.filtered')

    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    cmd = "java -Xmx12g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T MuTect2 -rf BadCigar -R {} {} {} " \
          " {} -o {} ".format(tmpdir, ref, coverage, args, input, vcf, cosmic)
    job_name = qsub.qsub('_'.join(('mutect',sample)), cmd, modules=modules, other='cpu=4|mem=16gb|walltime=120:00:00', waitlist=wait, **other_qsub_options)

    return job_name


if __name__ == '__main__':
    args = init()
    source = args.source
    outputdir = args.output
    log = args.log
    qsubdir = args.qsub
    dry = args.dry
    dbsnp = args.dbsnp
    coverage = args.coverage
    bed = args.bed
    cosmic = args.cosmic
    config = args.config

    if args.ini:
        pmgctools.read_vars(args.ini)

    if os.path.exists(config):
        with open(config, "r") as file:
            for line in file:
                sample = re.sub(' +', ',', line.rstrip()).replace(".bam", "").replace(".processed", "")
                if len(sample) > 0:
                    mutect(source=source, outputdir=outputdir, sample=sample, bed=bed, no_padding=args.no_padding,
                       dbsnp=dbsnp, cosmic=cosmic,
                       args=args.args, coverage=coverage, log=log, qsub=qsubdir, queue=args.queue, dry=dry)
    else:
        print("Error : couldn' find the config file : " + config)

