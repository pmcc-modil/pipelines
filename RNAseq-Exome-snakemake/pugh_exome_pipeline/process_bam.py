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
  Finds all the BAM in the source directory and run the recommended GATK pre-processing steps: Picard
  MarkDuplicates, realignment around indels, and base quality score recalibration. Options
  and filters chosen based on GATK best practices.
"""

import argparse
import os
import glob

import pmgctools
import qsub

KNOWN_1000G = '/mnt/work1/data/genomes/human/hg19/variantcallingdata/1000G_phase1.indels.hg19.vcf'
KNOWN_MILLS = '/mnt/work1/data/genomes/human/hg19/variantcallingdata/Mills_and_1000G_gold_standard.indels.hg19.vcf'


def init():
    parser = argparse.ArgumentParser(
        description='Finds all the BAM in the source directory and run the recommended GATK pre-processing steps')
    parser.add_argument('-s', '--source', required=True, help='Source directory')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-c', '--coverage',
                        help="coverage to downsample to at a given locus, enter an integer (e.g. 1000) or gatk to "
                             "use GATK default settings. If not provided, no downsampling will be performed.")
    parser.add_argument('-l', '--log', default='process.log', help='Log file name')
    parser.add_argument('-q', '--qsub', default="qsub", help='qsub directory')
    parser.add_argument('-D', '--dry-run', action='store_true', dest='dry', default=False,
                        help='dry run, will create qsub scripts but will not submit to the cluster')
    parser.add_argument('-I', '--ini', required=False, help='INI file.')
    parser.add_argument('-d', '--dbsnp', help='current dbSNP reference in VCF format')
    parser.add_argument('-g', '--bed', help='BED file containing intervals of interest')
    parser.add_argument('-S', '--species',
                        help='Species default=HUMAN. For non-human data extra Mills/1000Genomes indel resources, '
                             'which are provided by default to RealignerTargetCreator & BaseRecalibrator will not '
                             'be used.')
    parser.add_argument('-n', '--no-markduplicate', action='store_true', default=False, dest='no_markdup',
                        help='Do not run markduplicate.')
    parser.add_argument('-t', '--no-bqsr', action='store_true', default=False, dest='bqsr',
                        help='BQSR will not be run; needed for datasets < 150 Mbases')
    parser.add_argument('-C', '--config', default=None,
                        help='cocleaning config file. Put all samples for cocleaning in one line delimited with '
                             'space or ,')
    parser.add_argument('-Q', '--queue', default=None, help='Cluster queue you want to submit to.')
    options = parser.parse_args()
    return options


def mark_duplicate(source, outputdir, sample, waitlist=None, **other_qsub_options):
    """
    GATK Picard markduplicate
    :param source: source directory contains bam files need to be processed
    :param outputdir: output directory
    :param sample: sample name or bam file name to process
    :param waitlist: grid engine job waiting list
    :param other_qsub_options: other options for qsub script
    :return: job name, markdup_<sample>
    """
    tools = ['picard', 'samtools']
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    tmpdir = os.path.join(outputdir, pmgctools.tmpdir())

    if sample.endswith('.bam'):
        sample = sample[:-4]
    dedup = os.path.join(outputdir, sample + '.dedup.bam')
    metrics = os.path.join(outputdir, sample + '.dedup')
    picard_dir = 'picard_dir' if pmgctools.get_var('picard_dir') else 'EBROOTPICARD' #uhn or CC software

    cmd = 'mkdir -p {}\n'.format(tmpdir)
    cmd += "java -Xmx12g -Djava.io.tmpdir={} -jar ${}/picard.jar MarkDuplicates".format(tmpdir, picard_dir)
    cmd += " INPUT={} OUTPUT={} METRICS_FILE={} ASSUME_SORTED=true MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true".format(os.path.join(source, sample + '.bam'), dedup,
                                                                           metrics)
    cmd += "\nrm -rf {}".format(tmpdir)

    return qsub.qsub('markdup_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=4|mem=20gb|walltime=72:00:00', **other_qsub_options)


def indel_realignment(source, output, sample, bed=None, dbsnp=None, species='Human', coverage=None, smalldataset=False,
                      waitlist=None, extension='.dedup.bam', remove_source=True, **other_qsub_options):
    """
    GATK indel realignment
    :param source: source directory of the result of markduplicate (mark_duplicate result directory)
    :param sample: sample names to process, simple sample or multiple samples in list to co-clean
    :param bed: region file
    :param dbsnp: dbsnp vcf file
    :param species: HUMAN or not
    :param coverage: downsampling info, see source code
    :param smalldataset: small data set. do not run BSQR
    :param waitlist: job waiting list
    :param other_qsub_options: other options for qsub script
    :return: job name, Realignement_<sample>
    """
    tools = ['gatk', 'samtools']
    envs = pmgctools.check_vars(['REF'])
    if not bed:
        bed = pmgctools.get_var('bed')

    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    source = os.path.abspath(source)
    tmpdir = pmgctools.tmpdir()

    region = ""
    if bed:
        region = "--intervals " + os.path.abspath(bed) + " --interval_padding 100"
    ref = os.path.abspath(envs["REF"])
    if not dbsnp:
        dbsnp = os.path.abspath(pmgctools.check_vars(["dbSNP"])["dbSNP"])
    if species.upper() == "HUMAN":
        known_1000g = pmgctools.get_var("KNOWN_1000G")
        if known_1000g is None:
            known_1000g = KNOWN_1000G
        known_mills = pmgctools.get_var("KNOWN_MILLS")
        if known_mills is None:
            known_mills = KNOWN_MILLS
        known = "-known {} -known {}".format(known_1000g, known_mills)
    else:
        known = "-known {}".format(dbsnp)

    if not coverage:
        coverage = "-dt None"
    elif coverage is "gatk":
        coverage = ""
    elif coverage.isdigit():
        coverage = "-dcov " + coverage
    else:
        print("Error: wrong downsampling value " + coverage)

    input = "-I "
    if isinstance(sample, str):
        sample = [sample]

    input += ' -I '.join([os.path.join(source, i + extension) for i in sample])
    name = '_'.join(sample)[:200]
    intervals = name + '.intervals'

    # Start to work in the directory
    if not os.path.exists(output):
        os.makedirs(output)
    cmd = 'cd {}\nmkdir {}\n'.format(output, tmpdir)

    cmd += "java -Xmx8g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T RealignerTargetCreator --disable_auto_index_creation_and_locking_when_reading_rods -nt 4 -R {}" \
           " {} -o {} {} {} {}".format(tmpdir, ref, input, intervals, coverage, region, known)
    cmd += "\njava -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar --disable_auto_index_creation_and_locking_when_reading_rods -T IndelRealigner {} -nWayOut " \
           ".realigned.bam -targetIntervals {} -R {} {} {} -compress 0".format(tmpdir, input, intervals, ref, coverage,
                                                                               known)

    # remove intermediate files
    cmd += '\nrm -rf {}'.format(tmpdir)
    cmd += '\nrm {}'.format(intervals)
    new_extension = extension.replace('.bam', '.realigned.bam')
    for f in sample:
        if smalldataset:
            cmd += '\nmv {} {}'.format(f + new_extension, f + '.processed.bam')
            cmd += '\nmv {} {}'.format(f + new_extension.replace('.bam', '.bai'), f + '.processed.bai')
        if remove_source:
            cmd += '\nrm {}'.format(f + extension.replace('.bam', '.ba*'))

    cmd += '\n'

    return qsub.qsub('Realignment_' + name, cmd, modules=modules, waitlist=wait, other='cpu=8|mem=16gb|walltime=72:00:00', **other_qsub_options)


def bqsr(source, sample, bed=None, dbsnp=None, species='Human', coverage=None, extension='.dedup.realigned.bam',
         waitlist=None, **other_qsub_options):
    """
    GATK BQSR wrapper
    :param source: source directory of the result of indel_realignment (indel_realignment result directory)
    :param sample: sample names to process
    :param bed: region file
    :param dbsnp: :param dbsnp: dbsnp vcf file
    :param species: HUMAN or not
    :param coverage: downsampling info, see source code
    :param waitlist: job waiting list
    :param other_qsub_options: other options for qsub script
    :return: job name, BSQR_<sample>
    """
    tools = ['gatk']
    envs = pmgctools.check_vars(['REF'])
    modules = pmgctools.check_vars(tools)
    wait = ",".join(waitlist) if waitlist is not None else None
    tmpdir = os.path.join(source, pmgctools.tmpdir())

    if not bed:
        bed = pmgctools.get_var('bed')
    region = ""
    if bed:
        region = "--intervals " + bed + " --interval_padding 100"
    if not coverage:
        coverage = "-dt None"
    elif coverage is "gatk":
        coverage = ""
    elif coverage.isdigit():
        coverage = "-dcov " + coverage
    else:
        print("Error: wrong downsampling value " + coverage)
    recaldata = os.path.join(source, sample + ".recal_data.grp")
    ref = envs["REF"]
    if not dbsnp:
        dbsnp = pmgctools.check_vars(["dbSNP"])["dbSNP"]
    input = os.path.join(source, sample + extension)
    cmd = 'mkdir {}\n'.format(tmpdir)
    cmd += "java -Xmx4g -Djava.io.tmpdir={}  -jar $gatk_dir/GenomeAnalysisTK.jar -T BaseRecalibrator --disable_auto_index_creation_and_locking_when_reading_rods -nct 8 -I {} -o {}" \
          " -R {} -knownSites {} -rf BadCigar -cov ReadGroupCovariate -cov ContextCovariate -cov CycleCovariate -cov" \
          " QualityScoreCovariate {} {}".format(tmpdir, input, recaldata, ref, dbsnp, coverage, region)
#    if species.upper() == "HUMAN":
#        known_1000g = pmgctools.get_var("KNOWN_1000G")
#        if known_1000g is None:
#            known_1000g = KNOWN_1000G
#        known_mills = pmgctools.get_var("KNOWN_MILLS")
#        if known_mills is None:
#            known_mills = KNOWN_MILLS
#        cmd += " -knownSites {} -knownSites {}".format(known_1000g, known_mills)
    recal = os.path.join(source, sample + '.processed.bam')
    cmd += "\njava -Xmx4g -Djava.io.tmpdir={} -jar $gatk_dir/GenomeAnalysisTK.jar -T PrintReads --disable_auto_index_creation_and_locking_when_reading_rods -nct 8 -I {} -R {}" \
           " -BQSR {} -o {} -rf BadCigar {}".format(tmpdir, input, ref, recaldata, recal, coverage)
    cmd += '\nrm {} {} {}\n'.format(recaldata, input, input.replace('.bam', '.bai'))
    cmd += '\nrm -rf {}'.format(tmpdir)

    return qsub.qsub('BQSR_' + sample, cmd, modules=modules, waitlist=wait, other='cpu=8|mem=12gb|walltime=72:00:00', **other_qsub_options)


def process_bam(source, outputdir, samples, smalldataset=False, bed=None, dbsnp=None, species='Human',
                coverage=None, waitlist=None, no_markdup=False, **other_qsub_options):
    if isinstance(samples, str):
        samples = [samples]

    wait = None
    remove_source = False
    if not no_markdup:
        wait = []
        remove_source = True
        for sample in samples:
            wait.append(mark_duplicate(source, outputdir, sample, waitlist=waitlist, **other_qsub_options))
    else:
        wait = waitlist

    indel_wait = []
    dup_source = source if no_markdup else outputdir
    extension = '.bam' if no_markdup else '.dedup.bam'
    indel_wait.append(
        indel_realignment(dup_source, outputdir, samples, bed, dbsnp, species, coverage, smalldataset, waitlist=wait,
                          extension=extension, remove_source=remove_source, **other_qsub_options))
    bqsr_wait = []
    if not smalldataset:
        extension = '.dedup.realigned.bam' if not no_markdup else '.realigned.bam'
        for sample in samples:
            bqsr_wait.append(
                bqsr(outputdir, sample, bed, dbsnp, species, coverage, extension=extension, waitlist=indel_wait,
                     **other_qsub_options))

    return bqsr_wait + indel_wait


if __name__ == '__main__':
    args = init()
    if args.ini:
        pmgctools.read_vars(args.ini)
    pmgctools.check_vars(['gatk', 'samtools', 'picard', 'REF'])
    source, outputdir, config, log, qsubdir, dry, dbsnp, bed, species, coverage, smalldataset, config = args.source, \
                                                                                                        args.output, args.config, args.log, args.qsub, args.dry, args.dbsnp, args.bed, args.species, args.coverage, \
                                                                                                        args.bqsr, args.config
    # species command line > INI file, default: HUMAN
    if not species:
        species = pmgctools.get_var('SPECIES')
    if not species:
        species = 'HUMAN'

    remove_source = False if args.no_markdup else True

    if config is None:
        bams = glob.glob(os.path.join(source, '*.bam'))
        for bamfile in bams:
            sample = os.path.basename(bamfile)[:-4]  # remove .bam
            process_bam(source=source, outputdir=outputdir, samples=sample, smalldataset=smalldataset, bed=bed,
                        dbsnp=dbsnp, species=species, no_markdup=args.no_markdup,
                        coverage=coverage, log=log, qsub=qsubdir, dry=dry, queue=args.queue)

    else:
        with open(config) as f:
            for line in f:
                samples = line.rstrip().replace('.bam', '').replace('Sample_', '').replace(',', ' ').split()
                process_bam(source=source, outputdir=outputdir, samples=samples, smalldataset=smalldataset, bed=bed,
                            dbsnp=dbsnp, species=species, no_markdup=args.no_markdup,
                            coverage=coverage, log=log, qsub=qsubdir, dry=dry, queue=args.queue)
