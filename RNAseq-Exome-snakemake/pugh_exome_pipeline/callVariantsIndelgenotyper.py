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
This script runs IndelGenotyperV2 on tumor/normal pair of bam files
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
    parser.add_argument('-r', '--ref', default=None, help='Genome reference fasta file.')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-C', '--config', required=True,
                        help='mutect config file. Put tumor and normal sample pair in one line seperated by space or , '
                             'i.e tumor1.bam,normal1.bam')
    parser.add_argument('-a', '--args', default='', help='additional arguments to be provided to mutect')

    options = parser.parse_args()
    return options


def indelGenotyper(source, outputdir, sample, ref=None, args='', waitlist=None, extension=None, **other_qsub_options):
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())
    modules = pmgctools.check_vars(["java"])

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
    else:
        print("ERROR: can not find tomor/normal pair {}.".format(sample), file=sys.stderr)
        exit(1)
    if ref is None:
        ref = pmgctools.check_var("REF")

    sample = sample.replace(',', '__')[:200]
    vcf = os.path.join(outputdir, sample + '.indelgenotyper.vcf')
    metrics_file = os.path.join(outputdir, sample + '.indelgenotyper.outmetrics')

    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)

    cmd = "java -Xmx12g -Djava.io.tmpdir={} -jar {}/IndelGenotyper.36.3336-GenomeAnalysisTK.jar -T IndelGenotyperV2 -rf BadCigar --somatic --window_size 300 -R {} -o {} -metrics {} {} {}".format(tmpdir, os.path.dirname(__file__), ref, vcf, metrics_file, input, args)
    job_name = qsub.qsub('_'.join(['IndelGenotyper',sample]), cmd, modules=modules, other='cpu=1|mem=16gb|walltime=120:00:00', waitlist=waitlist, **other_qsub_options)

    return job_name


if __name__ == '__main__':
    args = init()

    if args.ini:
        pmgctools.read_vars(args.ini)

    if os.path.exists(args.config):
        with open(args.config, "r") as file:
            for line in file:
                sample = re.sub(' +', ',', line.rstrip()).replace(".bam", "").replace(".processed", "")
                if len(sample) > 0:
                    indelGenotyper(source=args.source, outputdir=args.output, sample=sample, ref=args.ref,
                       args=args.args, log=args.log, qsub=args.qsub, queue=args.queue, dry=args.dry)
    else:
        print("Error : couldn' find the config file : " + args.config)


