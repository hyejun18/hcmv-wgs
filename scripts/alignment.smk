import os

configfile: "scripts/config.yaml"


TARGETS = [
    # expand('data/processed/reads/{sample}.R1.flt.fq.gz', sample=config['samples']),
    # expand('data/processed/reads/{sample}.R2.flt.fq.gz', sample=config['samples']),
    # expand('data/processed/reads/{sample}.read_names.txt', sample=config['samples']),
    # expand('data/processed/reads/{sample}.R1.filtered.fq.gz', sample=config['samples']),
    # expand('data/processed/virus-stats/{sample}.jumps.txt.gz', sample=config['samples']['ID']),
    # expand('data/processed/alignments/{sample}.sorted.bam.bai', sample=config['samples']['ID']),
    expand('data/processed/alignments/{sample}.{species}.bam.bai', sample=config['samples']['ID'], species=['host', 'virus', 'discordant']),
]


rule all:
    input: TARGETS


rule align_to_genome:
    input: 
        fq1 = 'data/processed/reads/{sample}_1.trimmed.fq.gz',
        fq2 = 'data/processed/reads/{sample}_2.trimmed.fq.gz'
    output: 'data/processed/alignments/{sample}.bam'
    threads: 108
    params:
        genome_index = 'resources/refs/' + os.path.basename(config['reference_info']['reference_fasta'])[:-len('.fa.gz')] + \
            '.par_masked.ucsc.' + \
            config['reference_info']['viral_genome_assembly'],
    shell: 'bwa mem -t {threads} -M {params.genome_index} {input.fq1} {input.fq2} | samtools view -b -o {output} -'

rule sort_bam:
    input: 'data/processed/alignments/{sample}.bam'
    output: 'data/processed/alignments/{sample}.sorted.bam'
    threads: 12
    shell: 'samtools sort -@ {threads} -o {output} {input}'

rule index_bam:
    input: 'data/processed/alignments/{sample}.sorted.bam'
    output: 'data/processed/alignments/{sample}.sorted.bam.bai'
    threads: 12
    shell: 'samtools index -@ {threads} {input}'

rule collect_mapped_to_host_genome:
    input: 
        bam = 'data/processed/alignments/{sample}.sorted.bam',
        bam_idx = 'data/processed/alignments/{sample}.sorted.bam.bai'
    output: 'data/processed/alignments/{sample}.host.bam'
    threads: 1
    shell: 'samtools idxstats {input.bam} | cut -f 1 | \
            grep -v "GU937742.2" | \
            xargs samtools view -f 2 -F 256 -b -o {output} {input.bam}'

rule collect_mapped_to_viral_genome:
    input:
        bam = 'data/processed/alignments/{sample}.sorted.bam',
        bam_idx = 'data/processed/alignments/{sample}.sorted.bam.bai'
    output: 'data/processed/alignments/{sample}.virus.bam'
    threads: 1
    shell: 'samtools view -f 2 -F 256 -b -o {output} {input.bam} GU937742.2'

rule collect_discordant_reads:
    input: 
        bam = 'data/processed/alignments/{sample}.sorted.bam',
        bam_idx = 'data/processed/alignments/{sample}.sorted.bam.bai'
    output: 'data/processed/alignments/{sample}.discordant.bam'
    threads: 1
    shell: 'samtools view -F 258 -b -o {output} {input.bam}'

rule index_host_bam:
    input: 'data/processed/alignments/{sample}.host.bam'
    output: 'data/processed/alignments/{sample}.host.bam.bai'
    threads: 12
    shell: 'samtools index -@ {threads} {input}'

rule index_virus_bam:
    input: 'data/processed/alignments/{sample}.virus.bam'
    output: 'data/processed/alignments/{sample}.virus.bam.bai'
    threads: 12
    shell: 'samtools index -@ {threads} {input}'

rule index_discordant_bam:
    input: 'data/processed/alignments/{sample}.discordant.bam'
    output: 'data/processed/alignments/{sample}.discordant.bam.bai'
    threads: 12
    shell: 'samtools index -@ {threads} {input}'