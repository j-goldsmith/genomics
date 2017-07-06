import re
import urllib.request
import os.path
import subprocess

def sra_to_vcf(data_root, sra_name, ref):
    fastq = sra_to_fastq(sra_name, data_root)
    sam = fastq_to_sam(fastq, ref)
    bam = sam_to_bam(sam)
    #recalibrated_bam = recalibrate_bases(bam, ref)
    raw_variants = bam_to_raw_vcf(bam, ref)

    os.remove(fastq)
    os.remove(sam)
    os.remove(bam)

def bam_to_raw_vcf(bam,ref):
    filename = bam.replace('bam','vcf')

    if os.path.exists(filename) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/GenomeAnalysisTK.jar',
            '-T', 'HaplotypeCaller',
            '-R', ref,
            '-I', bam,
            '-o', filename])

    return filename

def sam_to_bam(sam):
    filename = sam.replace('sam','bam')
    sorted_f = filename[:-4] + '_sorted' + filename[-4:]
    rg_f = filename[:-4] + '_rg' + filename[-4:]

    if os.path.exists(sorted_f) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/picard.jar',
            'SortSam',
            'INPUT='+sam,
            'OUTPUT='+sorted_f,
            'SORT_ORDER=coordinate'])
        subprocess.run([
            'java',
            '-jar',
            '/bin/picard.jar',
            'AddOrReplaceReadGroups',
            'INPUT=' + sorted_f,
            'OUTPUT=' + rg_f,
            'RGID=1',
            'RGLB=libl',
            'RGPL=illumina',
            'RGSM='+sam,
            'RGPU=unit1'])

    if os.path.exists(filename) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/picard.jar',
            'MarkDuplicates',
            'INPUT=' + rg_f,
            'OUTPUT=' + filename,
            'CREATE_INDEX=true',
            'METRICS_FILE=metrics.txt'])

    os.remove(rg_f)
    os.remove(sorted_f)

    return filename

def sra_to_fastq(sra, data_root, gzip=True):
    in_file = os.path.join(data_root, 'sra', sra)
    out_dir = os.path.join(data_root,'fastq')
    out_file = os.path.join(out_dir, sra.replace('.sra','.fastq'))

    if os.path.exists(out_file) is False:
        subprocess.run([
            'fastq-dump',
            '-O',
            out_dir,
            in_file])

    return out_file

def fastq_to_sam(fastq, reference_genome, gzip=True):

    filename = fastq.replace('fastq','sam')
    if os.path.exists(filename) is False:
        subprocess.run([
            'bwa', 'mem',
            '-t', '6',
            '-M',
            reference_genome,
            fastq
           ], stdout=open(filename, 'w+'))

    return filename

def recalibrate_bases(bam, ref):
    table_f = bam[:-4] + '_recalibrated.table'
    post_table_f = bam[:-4] + '_post_recalibrated.table'
    filename = bam[:-4] + '_recalibrated.bam'

    if os.path.exists(table_f) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/GenomeAnalysisTK.jar',
            '-T', 'BaseRecalibrator',
            '-R', ref,
            '-I', bam,
            '-o', table_f])

    if os.path.exists(post_table_f) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/GenomeAnalysisTK.jar',
            '-T', 'BaseRecalibrator',
            '-R', ref,
            '-I', bam,
            '-BQSR', table_f,
            '-o', post_table_f])

        subprocess.run([
            'java',
            '-jar',
            '/bin/GenomeAnalysisTK.jar',
            '-T', 'PrintReads',
            '-R', ref,
            '-I', bam,
            '-BQSR', table_f,
            '-o', filename])

    return filename

if __name__ == "__main__":
    p = './phylos/data/sra/'
    for f in os.listdir(p):
        if os.path.isfile(os.path.join(p, f)):
            sra_to_vcf('./phylos/data/', f, "./phylos/data/ref/GCA_001865755.1_ASM186575v1_genomic.fna")
