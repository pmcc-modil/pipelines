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
  Finds all the VCF in the source directory and run the Ensembl variant_effect_predictor.
"""

import argparse
import os
import glob

import pmgctools
import qsub


def init():
    parser = argparse.ArgumentParser(
        description='Finds all the VCF in the source directory and run the Ensembl variant_effect_predictor.')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-e', '--extension', required=False, default='vcf', help='extension of the input files.')
    parser.add_argument('-o', '--output', required=False, help='Output directory')
    parser.add_argument('-S', '--species', required=False, default='homo_sapiens', help='species, default=home_sapiens')
    parser.add_argument('-A', '--assembly', required=True, help='Assembly, i.e. GRCh38.')
    parser.add_argument('-r', '--ref', required=False, help='Reference fasta file.')
    parser.add_argument('-f', '--format', required=False, default='vcf', choices=["ensembl", "vcf", "hgvs", "id"], help='Format of the input file, default=vcf.')
    parser.add_argument('-O', '--out-format', dest='out_format', required=False, default='vcf', choices=["vcf", "tab", "json"], help='Format of the output file, default=vcf.')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit to.')
    parser.add_argument('-a', '--args', default='', help='additional arguments to be provided to VEP')
    options = parser.parse_args()
    return options


def vep(vcf, outputdir, assembly, species='homo_sapiens', format='vcf', out_format='vcf', extension='vcf', ref=None, args=None, waitlist=None, **other_qsub_options):
    tools = ['vep']
    modules = pmgctools.check_vars(tools)

    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    sample = os.path.basename(vcf)[:-len(extension)-1]
    output = os.path.join(outputdir, sample+'.vep.'+extension)
    if not ref:
        ref = pmgctools.check_var('REF')

    cmd = 'variant_effect_predictor.pl --offline --quiet --force --no_stats --buffer_size 200 --fork 4 --species {} --assembly {} --dir_cache $cache_dir --input_file {} --output_file {} --format {} --fasta {} --{} {}'.format(species, assembly, vcf, output, format, ref, out_format, args)

    sub_dir = os.path.basename(os.path.dirname(vcf))
    if sub_dir == '.':
        sub_dir = ''
    return qsub.qsub('vep_' + sub_dir + '_' + sample, cmd, modules=modules, waitlist=waitlist, other='cpu=4|mem=16gb|walltime=72:00:00', **other_qsub_options)


if __name__ == '__main__':
    args = init()
    if args.ini:
        pmgctools.read_vars(args.ini)

    outputdir = args.output if args.output else args.source

    #vcfs = glob.glob(os.path.join(args.source, '*.vcf'))
    rel_start = len(args.source) if args.source.endswith('/') else len(args.source) + 1
    for dir_name, subdir_list, file_list in os.walk(args.source, followlinks=True):
        for file_name in file_list:
            if file_name.endswith(args.extension) and not file_name.endswith('vep.'+args.extension):
                if not os.path.exists(os.path.join(dir_name, file_name[:-len(args.extension)]+'vep.'+args.extension)):
                    vep(vcf=os.path.join(dir_name, file_name), outputdir=os.path.join(outputdir, dir_name[rel_start:]), species=args.species, assembly=args.assembly, format=args.format, out_format=args.out_format, extension=args.extension, ref=args.ref, args=args.args, log=args.log, qsub=args.qsub, dry=args.dry, queue=args.queue)

