configfile: 'scripts/config.yaml'

TARGETS = [
    # expand('data/processed/reads/{sample}.R1.flt.fq.gz', sample=config['samples']),
    # expand('data/processed/reads/{sample}.R2.flt.fq.gz', sample=config['samples']),
    # expand('data/processed/reads/{sample}.read_names.txt', sample=config['samples']),
    # expand('data/processed/reads/{sample}.R1.filtered.fq.gz', sample=config['samples']),
    # expand('data/processed/reads/{sample}.R2.filtered.fq.gz', sample=config['samples']),
    # expand('data/processed/annotations/{sample}.intersect.bed.gz', sample=config['samples']),
    # expand('data/processed/annotations/{sample}.intersect.bed.gz', sample=config['samples']),
    'data/processed/summary-stats/all.annotations.count.txt',
    # expand('data/processed/annotations/{sample}.unannotated.bam', sample=config['samples']['ID'])
]

rule all:
    input: TARGETS

#################### For unannotated BAM ####################
# rule bedtools_intersect_unannotated_R1:
#     input: 'data/processed/alignments/{sample}.host.sorted.bam'
#     output: 'data/processed/annotations/{sample}.R1.unannotated.bam'
#     threads: 2
#     params:
#         annoBed = 'resources/annotations/GRCh38.annotations.bed'
#     shell: 'samtools view -b -f 64 {input} | \
#             bedtools intersect -abam - -s -f 0.5 -b {params.annoBed} \
#             -v -split > {output}'

# rule bedtools_intersect_unannotated_R2:
#     input: 'data/processed/alignments/{sample}.host.sorted.bam'
#     output: 'data/processed/annotations/{sample}.R2.unannotated.bam'
#     threads: 2
#     params:
#         annoBed = 'resources/annotations/GRCh38.annotations.bed'
#     shell: 'samtools view -b -f 128 {input} | \
#             bedtools intersect -abam - -S -f 0.5 -b {params.annoBed} \
#             -v -split > {output}'

# rule merge_unannotated_bam:
#     input: r1 = 'data/processed/annotations/{sample}.R1.unannotated.bam',
#            r2 = 'data/processed/annotations/{sample}.R2.unannotated.bam'
#     output: 'data/processed/annotations/{sample}.unannotated.bam'
#     threads: 16
#     shell: 'samtools merge -@ {threads} {output} {input.r1} {input.r2}'
############################################################

#################### Annotation intersectBed ####################
rule bedtools_intersect:
    input: 'data/processed/alignments/{sample}.host.bam'
    output: 'data/processed/annotations/{sample}.intersect.bed.gz'
    threads: 1
    params:
        annoBed = 'resources/annotations/GRCh38.annotations.bed'
    shell: 'samtools view -b -F 256 {input} | \
            bedtools intersect -bed -abam - -f 0.5 -b {params.annoBed} \
            -wa -wb -split | gzip -c -> {output}'

rule summarize_annotations:
    input: 'data/processed/annotations/{sample}.intersect.bed.gz'
    output: 'data/processed/summary-stats/{sample}.annotated.count.txt'
    threads: 1
    shell: 'python3 scripts/summarize-annotations-pe.py {input} {output}'
############################################################

#################### Stats ####################
rule count_all_reads:
    input: 'data/reads/{sample}_1.fastq.gz'
    output: 'data/processed/summary-stats/{sample}.all.count.txt'
    threads: 1
    shell: "zcat {input} | awk 'NR % 4 == 1' | wc -l > {output}"

rule count_viral_reads:
    input: 'data/processed/alignments/{sample}.virus.bam'
    output: 'data/processed/summary-stats/{sample}.virus.count.txt'
    threads: 1
    shell: "samtools view -F 268 {input} | cut -f 1 | sort -T /hermes/hyejun/tmp | uniq | wc -l > {output}"

rule count_discordant_reads:
    input: 'data/processed/alignments/{sample}.discordant.bam'
    output: 'data/processed/summary-stats/{sample}.discordant.count.txt'
    threads: 1
    shell: "samtools view {input} | cut -f 1 | sort | uniq | wc -l > {output}"

rule count_host_reads:
    input: 'data/processed/alignments/{sample}.host.bam'
    output: 'data/processed/summary-stats/{sample}.host.count.txt'
    threads: 1
    shell: 'samtools view -F 268 {input} | cut -f 1 | sort -T /hermes/hyejun/tmp | uniq | wc -l > {output}'

rule count_trimmed_reads:
    input: 'data/processed/reads/{sample}_1.trimmed.fq.gz'
    output: 'data/processed/summary-stats/{sample}.trimmed.count.txt'
    threads: 1
    shell: "zcat {input} | awk 'NR % 4 == 1' | wc -l > {output}"

rule count_unmapped_reads:
    input: 'data/processed/alignments/{sample}.bam'
    output: 'data/processed/summary-stats/{sample}.unmapped.count.txt'
    threads: 1
    shell: "samtools view -f 12 {input} | cut -f 1 | uniq | wc -l > {output}"

rule merge_stat_counts:
    input: expand('data/processed/summary-stats/{sample}.virus.count.txt', sample=config['samples']['ID']),
           expand('data/processed/summary-stats/{sample}.unmapped.count.txt', sample=config['samples']['ID']),
           expand('data/processed/summary-stats/{sample}.discordant.count.txt', sample=config['samples']['ID']),
           expand('data/processed/summary-stats/{sample}.all.count.txt', sample=config['samples']['ID']),
           expand('data/processed/summary-stats/{sample}.trimmed.count.txt', sample=config['samples']['ID']),
           expand('data/processed/summary-stats/{sample}.host.count.txt', sample=config['samples']['ID']),
           expand('data/processed/summary-stats/{sample}.annotated.count.txt', sample=config['samples']['ID']),
    output: 'data/processed/summary-stats/all.annotations.count.txt'
    threads: 1
    shell: "python3 scripts/merge-annotations-count.py data/processed/summary-stats/ {output}"
############################################################
