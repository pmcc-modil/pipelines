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
Runs the GATK module HaplotypeCaller with the option to use the gVCF mode introduced in v3.0-0.
If --gvcf is specified, then callVariantsGenotypeGVCFS.py must be run to complete the variant
calling. Parameters for gVCF mode were chosen based on the [GATK workflow document]
(http://gatkforums.broadinstitute.org/discussion/3893/\
calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode). If --gvcf is
not specified, the variant calls are made as they were in GATK v2.x and are annotated with
GATK VariantAnnotator.

The --regions argument is STRONGLY recommended if running an exome or targeted experiement.
If --coverage is not specified GATK will downsample to a read depth of 250 for HaplotypeCaller
and 1000 for VariantAnnotataor, so it should be specified for high-coverage data.
"""

import os
import argparse
import sys
import re

import qsub
import pmgctools


def init():
    parser = argparse.ArgumentParser(
        description='Runs the GATK module HaplotypeCaller with the option to use the gVCF mode')
    parser.add_argument('-b', '--bam', required=True, help='bam file')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-r', '--reference', help=' FASTA reference sequence to which BAM was aligned')
    parser.add_argument('-d', '--dbsnp',
                        help='current dbSNP reference in VCF format (if sample on multiple lanes and need to merge)')
    parser.add_argument('-G', '--gvcf', action='store_true', default=False,
                        help='emits experimental reference confidence scores in GVCF mode for GenotypeGVCFs')
    parser.add_argument('-a', '--args',default='',
                        help='additional arguments to be provided to the HaplotypeCaller walker)'
                        )
    parser.add_argument('-g', '--bed',
                        help='BED file containing intervals of interest (if sample on multiple lanes and need to merge)'
                        )
    options = parser.parse_args()
    return options


def haplotype(bam, outputdir, ref=None, bed=None, dbsnp=None, gvcf=None, ext='.HaplotypeC.raw.vcf', args='', waitlist=None, **other_qsub_options):
    tools = ['gatk']
    modules = pmgctools.check_vars(tools)
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())

    wait = ",".join(waitlist) if waitlist is not None else None

    match = re.match('(.+?)(\.processed)*\.bam', os.path.basename(bam))
    if not match:
        print('ERROR: bam file name should end with ".bam" or ".processed.bam"', file=sys.stderr)
        exit(1)
    name = match.group(1)

    if not bed:
        bed = pmgctools.get_var('bed')
    if bed:
        args += " --intervals " + bed + " --interval_padding 100"

    if not ref:
        ref = pmgctools.check_var("REF")

    if not dbsnp:
        dbsnp = pmgctools.get_var("dbSNP")

    raw = os.path.join(outputdir, name + ext)

    cmd = 'mkdir -p {}\n'.format(tmpdir)
    cmd += "java -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller" \
          " {} -nct 4 -R {} -I {} -o {} -D {}  -stand_emit_conf 10 -rf BadCigar ".\
        format(tmpdir, args, ref, bam, raw, dbsnp)
    if gvcf:
        cmd+=" -ERC GVCF -variant_index_type LINEAR -variant_index_parameter 128000 "
    cmd += '\nrm -rf {}'.format(tmpdir)

    return qsub.qsub('haplotypeCaller_' + name, cmd, modules=modules, other='cpu=4|mem=8gb|walltime=72:00:00', **other_qsub_options)


if __name__ == '__main__':
    arguments = init()
    bam = arguments.bam
    outputdir = arguments.output
    log = arguments.log
    qsubdir = arguments.qsub
    dry = arguments.dry
    dbsnp = arguments.dbsnp
    ref = arguments.reference
    bed = arguments.bed
    gvcf = arguments.gvcf
    args=arguments.args

    haplotype(bam=bam, outputdir=outputdir, gvcf=gvcf, ref=ref, bed=bed, dbsnp=dbsnp, args=args, log=log,
               qsub=qsubdir, dry=dry)
