import re
import urllib.request
import os.path
import subprocess

def sra_to_vcf(data_root, sra_name, ref):
    fastq = sra_to_fastq(sra_name, data_root)
    sam = fastq_to_sam(fastq, ref)
    bam = sam_to_bam(sam)
    recalibrated_bam = recalibrate_bases(bam, ref)

def sam_to_bam(sam):
    filename = sam.replace('sam','bam')
    sorted_f = filename[:-4] + '_sorted' + filename[-4:]

    if os.path.exists(sorted_f) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/picard.jar',
            'SortSam',
            'INPUT='+sam,
            'OUTPUT='+sorted_f,
            'SORT_ORDER=coordinate'])

    if os.path.exists(filename) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/picard.jar',
            'MarkDuplicates',
            'INPUT=' + sorted_f,
            'OUTPUT=' + filename,
            'METRICS_FILE=metrics.txt'])

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
            '-L', '20',
            '-o', table_f])

    if os.path.exists(post_table_f) is False:
        subprocess.run([
            'java',
            '-jar',
            '/bin/GenomeAnalysisTK.jar',
            '-T', 'BaseRecalibrator',
            '-R', ref,
            '-I', bam,
            '-L', '20',
            '-BQSR', table_f,
            '-o', post_table_f])

        subprocess.run([
            'java',
            '-jar',
            '/bin/GenomeAnalysisTK.jar',
            '-T', 'PrintReads',
            '-R', ref,
            '-I', bam,
            '-L', '20',
            '-BQSR', table_f,
            '-o', filename])

    return filename

if __name__ == "__main__":
    sra_to_vcf('./phylos/data/', 'SRR4448725.sra', "./phylos/data/ref/GCA_001865755.1_ASM186575v1_genomic.fna")
