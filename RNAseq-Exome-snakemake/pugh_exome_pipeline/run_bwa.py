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
This script runs BWA on gzipped FASTQ files, which are found by searching the source directory for files ending in "*R1.fastq.gz"
"""

import re
import sys
import glob
import os
import argparse

import qsub
import pmgctools
import process_bam


def init():
    parser = argparse.ArgumentParser(
        description='This script runs BWA on gzipped FASTQ files, which are found by searching the source '
                    'directory for files ending in "*R1.fastq.gz"')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-p', '--project', required=True, help='Project name')
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsubdir', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-d', '--dbsnp',
                        help='current dbSNP reference in VCF format (if sample on multiple lanes and need to merge)')
    parser.add_argument('-g', '--bed',
                        help='BED file containing intervals of interest (if sample on multiple lanes and need to merge)')
    parser.add_argument('-n', '--no-markduplicate', action='store_true', default=False, dest='no_markdup',
                        help='Do not run markduplicate.')
    parser.add_argument('-S', '--species', default="HUMAN",
                        help='Species (if sample on multiple lanes and need to merge)')
    parser.add_argument('-r', '--bwa-index', dest='bwaindex', default=None, help='Indexed reference genome')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit to.')
    parser.add_argument('-L', '--no-lane-process', action='store_false', dest="lane_process", default=True,
                        help=' will not process lanes, will just merge bam files')
    parser.add_argument('-O', '--old-format', action='store_true', dest="old_format", default=False,
                        help='Each sample is in a separate directory.')
    parser.add_argument('-C', '--sample-list', dest="sample_list", help='Sample list file. Each line list one sample name. This is only necessary if samples are on multiple lanes.')
    args = parser.parse_args()
    return args


def bwa(sample, output, read1, read2=None, index=None, readgroup=None, **other_qsub_options):
    """
    generate script for bwa and qsub to the cluster.
    :param sample: sample name
    :param output: output file prefix with directory, like /home/myname/bam/my_sample
    :param read1: R1 file name
    :param read2: R2 file name if pair end
    :param index: bwa index prefix, if None get from environment variable BWAINDEX
    :param readgroup: Complete read group header line
    :param other_qsub_options: other options to pass to qsub. please refer to pmgctools.qsub
    :return: qsub job name, "bwa_<sample>"
    """

    if read2 is None: read2 = ''
    tools = ['bwa', 'samtools']
    modules = pmgctools.check_vars(tools)

    if index is None:
        pmgctools.check_vars(['BWAINDEX'])
        index = pmgctools.VARS['BWAINDEX']

    cmd = 'bwa mem -M -t4'
    if readgroup is not None:
        cmd += ' -R ' + readgroup
    cmd += " {} {} {}>{}.sam\n".format(index, read1, read2, output)
    #cmd += "samtools view -bhS {0}.sam | samtools sort -@4 - {0}\n".format(output)
    cmd += "samtools sort -@4 -O bam -o {0}.bam -T {0} {0}.sam\n".format(output)
    cmd += "samtools index {}.bam\n".format(output)
    cmd += "rm {}.sam\n".format(output)

    name = 'bwa_' + sample
    return qsub.qsub(name, cmd, modules=modules, other='cpu=4|mem=8gb|walltime=72:00:00', **other_qsub_options)


def merge_bam(datadir, sample, samplelist, extension=".bam", waitlist=None, **other_qsub_options):
    """
    merge bam files if sample on multiple lanes
    :param datadir: the directory contains bam files to merge
    :param sample: sample name/bam file name after merge
    :param samplelist: bam file names to merge
    :param lane_process: if the bam files were processed
    :param waitlist: job waiting list
    :param other_qsub_options: other oprions for qsub script
    :return: job name MergeSam_<sample>
    """
    tools = ['picard', 'samtools']
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist)
    tmpdir = os.path.join(datadir, pmgctools.tmpdir())
    output = os.path.join(datadir, sample + '.bam')
    picard_dir = 'picard_dir' if pmgctools.get_var('picard_dir') else 'EBROOTPICARD' #uhn or CC software

    cmd = 'mkdir -p {}\n'.format(tmpdir)
    cmd += "java -Xmx8g -Djava.io.tmpdir={} -jar ${}/picard.jar MergeSamFiles  ".format(tmpdir, picard_dir)

    samples = [os.path.join(datadir, i + extension) for i in samplelist]

    cmd += "INPUT={} OUTPUT={} ASSUME_SORTED=true USE_THREADING=true".format(' INPUT='.join(samples), output)
    cmd += "\nsamtools index {}\n".format(output)
    #remove the cmobined files, including .bam, .dedup, .dedup.realigned.bam
    cmd += "rm -rf {}\n".format(tmpdir)
    cmd += "cd {}\nrm {}.*\n".format(datadir, '.* '.join(samplelist))

    return qsub.qsub('MergeSam_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=4|mem=12gb|walltime=24:00:00', **other_qsub_options)


def bwa_sample(source, outputdir, project, samplename=None, index=None, lane_process=True, no_markdup=False, bed=None,
               dbsnp=None, species='HUMAN', **other_qsub_options):
    sample_id = samplename

    waitlist = []

    # create output dir if not exists
    if not os.path.exists(outputdir):
        os.makedirs(outputdir, 0o770)

    fastqlist = glob.glob(os.path.join(source, sample_id + "_*R1.fastq*")) + glob.glob(os.path.join(source, sample_id + "_*R1_001.fastq*"))
    multiple = True if len(fastqlist) > 1 else False
    samplelist = []
    extension = ".bam" #for merge_bam function
    for read1 in fastqlist:
        pattern = re.compile('.*/' + sample_id + '_([A-Z0-9]+)_([LS]\d+)_R1(_001)*.fastq*')
        matched = pattern.match(read1)
        lane, barcode = '', ''
        if matched:
            barcode = matched.group(1)
            lane = matched.group(2)

        if 'R1_001.fastq' in read1:
            read2 = read1.replace('R1_001.fastq', 'R2_001.fastq')
        else:
            read2 = read1.replace('R1.fastq', 'R2.fastq')
        if not os.path.exists(read2):
            read2 = None

        if multiple:
            sample = os.path.basename(read1).split('_R1.fastq')[0]
        else:
            sample = sample_id

        samplelist.append(sample)
        readgroup = r"'@RG\tID:1\tSM:{}\tPL:Illumina\tPU:{}.{}\tLB:{}'".format(sample_id, barcode, lane, project)

        job_name = bwa(sample=sample, output=os.path.join(outputdir, sample), read1=read1, read2=read2, index=index,
                       readgroup=readgroup, **other_qsub_options)
        waitlist.append(job_name)

        #lane level process bam
        if multiple and lane_process:
            if no_markdup:
                extension = '.realigned.bam'
                indel_job = process_bam.indel_realignment(source=outputdir, output=outputdir, sample=sample, bed=bed, dbsnp=dbsnp,
                                                      species=species, coverage=None, smalldataset=False, extension='.bam',
                                                      waitlist=[job_name], **other_qsub_options)
            else:
                extension = '.dedup.realigned.bam'
                markdup_job = process_bam.mark_duplicate(source=outputdir, outputdir=outputdir, sample=sample,
                                                     waitlist=[job_name], **other_qsub_options)
                indel_job = process_bam.indel_realignment(source=outputdir, output=outputdir, sample=sample, bed=bed, dbsnp=dbsnp,
                                                      species=species, coverage=None, smalldataset=False,
                                                      waitlist=[markdup_job], **other_qsub_options)
            waitlist.append(indel_job)

    #merge bam if sample on multiple lanes
    if multiple:
        job_name = merge_bam(datadir=outputdir, sample=sample_id, samplelist=samplelist, extension=extension, waitlist=waitlist,
                             **other_qsub_options)

    return job_name


if __name__ == '__main__':
    args = init()
    if args.ini:
        pmgctools.read_vars(args.ini)
    pmgctools.check_vars(['bwa', 'samtools', 'BWAINDEX'])
    if args.old_format:
        for sample in glob.glob(os.path.join(args.source, 'Sample_*')):
            if os.path.isdir(sample):
                bwa_sample(source=sample, outputdir=args.output, project=args.project, index=args.bwaindex, log=args.log,
                           qsub=args.qsubdir, dry=args.dry, queue=args.queue, lane_process=args.lane_process, no_markdup=args.no_markdup,
                           bed=args.bed, dbsnp=args.dbsnp, species=args.species)
    else:
        if args.sample_list:
            with open(args.sample_list) as f:
                for line in f:
                    bwa_sample(source=args.source, outputdir=args.output, project=args.project, samplename=line.strip(), index=args.bwaindex, log=args.log,
                           qsub=args.qsubdir, dry=args.dry, queue=args.queue, lane_process=args.lane_process, no_markdup=args.no_markdup,
                           bed=args.bed, dbsnp=args.dbsnp, species=args.species)
        else:
            pattern = re.compile('.*/(.+)_([A-Z0-9]+)_([LS]\d+)_R1(_001)*.fastq*')
            sample_list = set()
            for fastq in glob.glob(os.path.join(args.source, "*R1.fastq*")) + glob.glob(os.path.join(args.source, "*R1_001.fastq*")):
                matched = pattern.match(fastq)
                if matched:
                    sample = matched.group(1)
                else:
                    splitor = '_R1_001.fastq' if '_R1_001.fastq' in fastq else '_R1.fastq'
                    sample = os.path.basename(fastq).split(splitor)[0]
                sample_list.add(sample)

            for sample in sample_list:
                bwa_sample(source=args.source, outputdir=args.output, project=args.project, samplename=sample, index=args.bwaindex, log=args.log,
                           qsub=args.qsubdir, dry=args.dry, queue=args.queue, lane_process=args.lane_process, no_markdup=args.no_markdup,
                           bed=args.bed, dbsnp=args.dbsnp, species=args.species)

