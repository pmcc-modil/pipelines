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
This script is to merge gvcf files.
"""

import glob
import os
import argparse
import sys

import pmgctools
import qsub

def init():
    parser = argparse.ArgumentParser(description='Merge gvcf files.')
    parser.add_argument('-s', '--source', required=True, help='source directories for vcf files')
    parser.add_argument('-i', '--input', nargs='*', help='gvcf files to merge')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-n', '--name', help='Output vcf name')
    parser.add_argument('-C', '--conf', help='conf file, each line contains result file name=gvc1f.vcf+gvcf2.vcf+...')
    parser.add_argument('-l', '--log', default='process.log', help='Log file.')
    parser.add_argument('-q', '--qsub', default='qsub', help='Directory to store the qusb script.')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False, help='Dry run.')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit to.')

    return parser.parse_args()


def merge_gvcf(input: list, output: str, **other_qsub_options) -> str:
    tools = ['gatk']
    dbsnp = pmgctools.get_var('dbSNP')
    ref = pmgctools.check_var('REF')
    modules = pmgctools.check_vars(tools)
    tmpdir = os.path.join(os.path.dirname(output), pmgctools.tmpdir())
    name = os.path.splitext(os.path.basename(output))[0]

    cmd = 'mkdir -p {}\n'.format(tmpdir)
    cmd += 'export _JAVA_OPTIONS="-XX:ParallelGCThreads=2"\n'
    cmd += "java -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T GenotypeGVCFs -R {} -o {} ".format(tmpdir, ref, output)
    if dbsnp:
        cmd += '-D {} '.format(dbsnp)
    cmd += '--variant ' + ' --variant '.join(input)

    cmd += "\nrm -rf {}\n".format(tmpdir)

    return qsub.qsub('GenotypeGVCFs_' + name, cmd, modules=modules, other='cpu=1|mem=8gb|walltime=24:00:00', **other_qsub_options)



if __name__ == '__main__':
    args = init()

    source = args.source if args.source else ''

    if args.conf:
        with open(args.conf, 'r') as f:
            for line in f:
                output_file = line.split('=')[0]
                gvcfs = line.strip().split('=')[1].split('+')
                gvcfs = [os.path.join(source, i) for i in gvcfs]
                merge_gvcf(input=gvcfs, output=os.path.join(args.output, output_file), qsub=args.qsub, dry=args.dry, queue=args.queue, log=args.log)
    else:
        if not args.name:
            raise NameError('result vcf name not defined.')

        gvcfs = [os.path.join(source, i) for i in args.input]
        output_file = args.name if args.name.endswith('.vcf') else args.name + '.vcf'

        merge_gvcf(input=gvcfs, output=os.path.join(args.output, output_file),
                   qsub=args.qsub, dry=args.dry, queue=args.queue, log=args.log)

