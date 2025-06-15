import os
# import subprocess as sp

configfile: "scripts/config.yaml"


MAX_THREADS = os.cpu_count()


TARGETS = [
    expand('data/processed/reads/{sample}_{read}trimmed.fq.gz', sample=config['samples']['ID'], read=['1.', '2.']),
]


rule all:
    input: TARGETS


rule overall_preprocessing_without_merge:
    input: 
        fq1 = 'data/reads/{sample}_1.fastq.gz',
        fq2 = 'data/reads/{sample}_2.fastq.gz'
    output: 
        fq1 = 'data/processed/reads/{sample}_1.trimmed.fq.gz',
        fq2 = 'data/processed/reads/{sample}_2.trimmed.fq.gz',
    threads: 16
    log: 'logs/{sample}.fastp.html'
    shell: 'fastp \
                -i {input.fq1} -I {input.fq2} \
                -o {output.fq1} -O {output.fq2} \
                --thread {threads} \
                -h {log} -j /dev/null \
                --thread {threads}'
