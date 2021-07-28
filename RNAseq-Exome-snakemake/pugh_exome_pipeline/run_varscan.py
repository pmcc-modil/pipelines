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
This script runs varscan on tumor/normal pair of bam files
"""

import sys
import os
import argparse
import re

import qsub
import pmgctools


def init():
    parser = argparse.ArgumentParser(
        description='This script runs varscan on tumor/normal pair of bam files')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-g', '--bed', help='BED file containing intervals of interest')
    parser.add_argument('-Q', '--queue', help='grid engine queue to use')
    parser.add_argument('-t', '--tool', choices=['somatic', 'mpileup2cns', 'copynumber'], help='sub-command for varscan')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-v', '--vcf', required=False, action='store_true', help='output in vcf format.')
    parser.add_argument('-a', '--args',default='',
                        help='additional arguments to be provided to varscan.)'
                        )
    parser.add_argument('-C', '--config', required=True,
                        help='varscan config file. Put tumor and normal sample pair in one line seperated by space or , '
                             'i.e tumor1.bam,normal1.bam')

    options = parser.parse_args()
    return options


def somatic(source, outputdir, tumor, normal, bed=None, args='', vcf=False, waitlist=None, extension=None, **other_qsub_options):
    tools = ['samtools', 'varscan']
    modules = pmgctools.check_vars(tools)
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())

    if not extension:
        if os.path.exists(os.path.join(source, tumor + '.processed.bam')):
            extension = '.processed.bam'
        elif os.path.exists(os.path.join(source, tumor + '.bam')):
            extension = '.bam'
        else:
            print('ERROR: I can not find bam files for samples "{}"'.format(tumor + ',' + normal), file=sys.stderr)
            exit(1)

    ref = pmgctools.check_var("REF")

    sample = tumor + '__' + normal
    sample = sample[:200]
    output = os.path.join(outputdir, sample)

    cmd = 'mkdir -p {0}\nexport TMPDIR={0}\n'.format(tmpdir)
    cmd += 'set -o pipefail\n'
    cmd += "samtools mpileup -B -q1 -d10000000"
    if not bed:
        bed = pmgctools.get_var('bed')
    if bed:
        cmd += ' -l {}'.format(bed)

    normal_bam = os.path.join(source, normal + extension)
    tumor_bam = os.path.join(source, tumor + extension)
    cmd += ''' -f {} {} {} | awk -F"\\t" '$4 > 0 && $7 >0' | '''.format(ref, normal_bam, tumor_bam)
    cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $varscan_dir/VarScan.jar somatic - {} --mpileup 1 {}".format(
        tmpdir, output, args)
    if vcf:
        cmd += ' --output-vcf 1'
    cmd += '\nrm -rf {}'.format(tmpdir)

    return qsub.qsub('varscan_' + sample, cmd, modules=modules, other='cpu=3|mem=20gb|walltime=72:00:00', waitlist=waitlist,
                     **other_qsub_options)


def single(source, outputdir, sample, bed=None, vcf=False, waitlist=None, extension=None, **other_qsub_options):
    tools = ['samtools', 'varscan', 'python3']
    modules = pmgctools.check_vars(tools)
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())
    wait = ",".join(waitlist) if waitlist is not None else None

    if not extension:
        if os.path.exists(os.path.join(source, sample + '.processed.bam')):
            extension = '.processed.bam'
        elif os.path.exists(os.path.join(source, sample + '.bam')):
            extension = '.bam'
        else:
            print('ERROR: I can not find bam files for sample "{}"'.format(sample), file=sys.stderr)
            exit(1)

    ref = pmgctools.check_var("REF")

    output = os.path.join(outputdir, sample)

    cmd = 'mkdir -p {0}\nexport TMPDIR={0}\n'.format(tmpdir)
    cmd += "samtools mpileup -B -q1 -d10000000"
    if not bed:
        bed = pmgctools.get_var('bed')
    if bed:
        cmd += ' -l {}'.format(bed)

    sample_bam = os.path.join(source, sample + extension)
    cmd += ''' -f {} {} | '''.format(ref, sample_bam, output)

    if vcf:
        cmd += "java -Xmx8g -Djava.io.tmpdir={} -jar $varscan_dir/VarScan.jar mpileup2cns - --output-vcf 1 > {}.cns.vcf".format(
            tmpdir, output, output)
    else:
        cmd += "java -Xmx8g -Djava.io.tmpdir={} -jar $varscan_dir/VarScan.jar mpileup2cns - > {}.cns".format(
            tmpdir, output, output)

    return qsub.qsub('varscan_' + sample, cmd, modules=modules, other='cpu=4|mem=12gb|walltime=24:00:00', waitlist=wait,
                     **other_qsub_options)


def copynumber(source, outputdir, tumor, normal, bed=None, args='', waitlist=None, extension=None, **other_qsub_options):
    tools = ['samtools', 'varscan']
    modules = pmgctools.check_vars(tools)
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())

    if not extension:
        if os.path.exists(os.path.join(source, tumor + '.processed.bam')):
            extension = '.processed.bam'
        elif os.path.exists(os.path.join(source, tumor + '.bam')):
            extension = '.bam'
        else:
            print('ERROR: I can not find bam files for samples "{}"'.format(tumor + ',' + normal), file=sys.stderr)
            exit(1)

    ref = pmgctools.check_var("REF")

    sample = tumor + '__' + normal
    sample = sample[:200]
    output = os.path.join(outputdir, sample)

    cmd = 'mkdir -p {0}\nexport TMPDIR={0}\n'.format(tmpdir)
    cmd += 'set -o pipefail\n'
    cmd += "samtools mpileup -B -q1 -d10000000"
    if not bed:
        bed = pmgctools.get_var('bed')
    if bed:
        cmd += ' -l {}'.format(bed)

    normal_bam = os.path.join(source, normal + extension)
    tumor_bam = os.path.join(source, tumor + extension)
    cmd += ''' -f {} {} {} | awk -F"\\t" '$4 > 0 && $7 >0' | '''.format(ref, normal_bam, tumor_bam, output)
    cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar $varscan_dir/VarScan.jar copynumber - {} --mpileup 1 {}\n".format(
        tmpdir, output, args)
    cmd += '\nrm -rf {}'.format(tmpdir)

    return qsub.qsub('varscan_' + sample, cmd, modules=modules, other='cpu=3|mem=20gb|walltime=72:00:00', waitlist=waitlist,
                     **other_qsub_options)


if __name__ == '__main__':
    args = init()

    if args.ini:
        pmgctools.read_vars(args.ini)

    if os.path.exists(args.config):
        with open(args.config, "r") as file:
            for line in file:
                sample = re.sub(' +', ',', line.strip()).replace(".bam", "")
                samples = sample.split(',')
                if not args.tool or args.tool in ['somatic', 'mpileup2cns']:
                    if len(samples) == 2:
                        tumor, normal = samples
                        somatic(source=args.source, outputdir=args.output, tumor=tumor, normal=normal, args=args.args, vcf=args.vcf,
                                bed=args.bed, log=args.log, qsub=args.qsub, dry=args.dry, queue=args.queue)
                    elif len(samples) == 1:
                        single(source=args.source, outputdir=args.output, sample=samples[0], bed=args.bed, vcf=args.vcf, log=args.log, qsub=args.qsub, dry=args.dry, queue=args.queue)
                    else:
                        print("Error: I can not find sample or pair for '{}'".format(line.strip()), file=sys.stderr)
                elif args.tool == 'copynumber':
                    if len(samples) == 2:
                        tumor, normal = samples
                        copynumber(source=args.source, outputdir=args.output, tumor=tumor, normal=normal,args=args.args,
                                bed=args.bed, log=args.log, qsub=args.qsub, dry=args.dry, queue=args.queue)

    else:
        print("Error : couldn' find the config file : " + args.config, file=sys.stderr)
