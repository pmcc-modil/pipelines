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
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-c', '--strelka-config', dest='stralka_config', required=False, help='Stralka user configuration file.')
    parser.add_argument('-r', '--ref', required=False, help='genome reference fasta file.')
    parser.add_argument('-C', '--config', required=True,
                        help='sample config file. Put tumor and normal sample pair in one line seperated by space or , '
                             'i.e tumor1.bam,normal1.bam')

    options = parser.parse_args()
    return options


def strelka(outputdir, tumor, normal, ref=None, config=None, waitlist=None, **other_qsub_options):

    tools = ['strelka']
    modules = pmgctools.check_vars(tools)

    if not ref:
        ref = pmgctools.check_var('REF')

    if not config:
        config = os.path.join(os.path.dirname(__file__), 'strelka_config_bwa.ini')

    sample = os.path.basename(tumor).replace(".bam", "").replace(".processed", "") + '__' + os.path.basename(normal).replace(".bam", "").replace(".processed", "")
    sample = sample[:200]

    result_dir = os.path.join(outputdir, sample)
    tmpdir = os.path.join(result_dir, pmgctools.tmpdir())

    if not os.path.exists(result_dir):
        os.makedirs(result_dir, 0o770)

    cmd = 'configureStrelkaWorkflow.pl --normal={} --tumor={} --ref={} --config={} --output-dir={}\n'.format(normal, tumor, ref, config, tmpdir)
    cmd += 'make -j 4 -C {}\n'.format(tmpdir)
    cmd += 'mv {}/* {}\n'.format(os.path.join(tmpdir, 'results'), result_dir)
    cmd += 'rm -rf {}\n'.format(tmpdir)


    return qsub.qsub('strelka_' + sample, cmd, modules=modules, other='cpu=4|mem=8gb|walltime=72:00:00', waitlist=waitlist, **other_qsub_options)


if __name__ == '__main__':
    args = init()

    if args.ini:
        pmgctools.read_vars(args.ini)

    if os.path.exists(args.config):
        with open(args.config, "r") as file:
            for line in file:
                sample = re.sub(' +', ',', line.rstrip()).replace(".bam", "").replace(".processed", "")
                pair = sample.split(",")
                if os.path.exists(os.path.join(args.source, pair[0] + '.processed.bam')):
                    extension = '.processed.bam'
                else:
                    extension = '.bam'

                tumor = os.path.join(args.source, pair[0] + extension)
                normal = os.path.join(args.source, pair[1] + extension)

                if not (os.path.exists(tumor) and os.path.exists(normal)):
                    print('ERROR: I can not find bam files for samples "{}"'.format(sample), file=sys.stderr)

                strelka(outputdir=args.output, tumor=tumor, normal=normal, config=args.stralka_config, ref=args.ref, log=args.log, qsub=args.qsub, dry=args.dry)
    else:
        print("Error : couldn' find the config file : " + args.config)
