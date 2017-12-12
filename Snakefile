shell.prefix("set -eo pipefail; ")

configfile: "config.yaml"
localrules: all

import glob


SAMPLES = glob.glob('samples/*.fastq.gz')
SAMPLES = [sample.replace('.fastq.gz','') for sample in SAMPLES]
SAMPLES = [sample.replace('samples/','') for sample in SAMPLES]


ALL_FASTQ = expand('samples/{sample}.fastq.gz', sample=SAMPLES)
ALL_TRIMMED_FASTQ = expand('trimmed_reads/{sample}.trimmed.fastq.gz', sample=SAMPLES)

rule all:
    input: 
        ALL_TRIMMED_FASTQ,
        expand('qc/{sample}_fastqc.html', sample=SAMPLES),
        expand('qc/{sample}_fastqc.zip', sample=SAMPLES),
        expand('qc/{sample}.trimmed_fastqc.html', sample=SAMPLES),
        expand('qc/{sample}.trimmed_fastqc.zip', sample=SAMPLES),
        expand('mapped/bams/{sample}.bam', sample=SAMPLES),
        expand('mapped/bams/{sample}.sortedByName.bam', sample=SAMPLES),
        expand('mapped/counts_strict/{sample}.counts.tsv', sample=SAMPLES),
        expand('mapped/counts_strict/{sample}.counts.noversion.tsv', sample=SAMPLES),
        expand('mapped/featureCounts/fcounts.tsv', sample=SAMPLES),
        'mapped/featureCounts/fcounts.noversion.tsv',
        expand('rsem/{sample}.isoforms.results', sample=SAMPLES),
        expand('rsem/{sample}.genes.results', sample=SAMPLES),
        expand('rsem/{sample}.genome.bam', sample=SAMPLES),
        expand('rsem/{sample}.genome.sorted.bam', sample=SAMPLES),
        expand('rsem/{sample}.genome.sorted.bam.bai', sample=SAMPLES)  
        

rule qc_raw:
    input: 'samples/{sample}.fastq.gz'
    params:
        out_dir = 'qc'
    output:
       'qc/{sample}_fastqc.html',
       'qc/{sample}_fastqc.zip'
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input}
        '''
     
rule trim_se:
     input: "samples/{sample}.fastq.gz"
     output: "trimmed_reads/{sample}.trimmed.fastq.gz"
     threads: 4
     shell: 
         """
         java -jar {config[path]}trimmomatic-0.36.jar SE -phred33 -threads {threads} \
         {input} {output} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
         """
rule qc_trimmed:
    input: 'trimmed_reads/{sample}.trimmed.fastq.gz'
    params:
        out_dir = 'qc'
    output:
       'qc/{sample}.trimmed_fastqc.html',
       'qc/{sample}.trimmed_fastqc.zip'
    shell:
        r'''
            fastqc -o {params.out_dir} -f fastq {input}
        '''

rule rsem:
    input: 'trimmed_reads/{sample}.trimmed.fastq.gz'
    output: 
        'rsem/{sample}.isoforms.results',
        'rsem/{sample}.genes.results',
        'rsem/{sample}.genome.bam',
        'rsem/{sample}.genome.sorted.bam',
        'rsem/{sample}.genome.sorted.bam.bai'              
        
    threads: 24
    shell:
        """
        rsem-calculate-expression --star --star-path {config[star]} --star-gzipped-read-file --output-genome-bam --sort-bam-by-coordinate -p {threads} {input} {config[rsem]} rsem/{wildcards.sample}
        """
        
rule map_star:
    input: 'trimmed_reads/{sample}.trimmed.fastq.gz'
    output: 'mapped/bams/{sample}.bam'
    params:
        prefix = 'mapped/bams/{sample}',
        unmapped = 'unmapped/fastq/{sample}',
        starlogs = 'mapped/starlogs'
    threads: 16
    shell:
        '''
        STAR --runThreadN {threads} --genomeDir {config[index]} --outFileNamePrefix {params.prefix} --readFilesIn {input} --outSAMtype BAM SortedByCoordinate --outFilterMatchNmin 50 --outFilterMismatchNmax 100 --readFilesCommand zcat --outReadsUnmapped {params.unmapped}  
        mv {params.prefix}Aligned.sortedByCoord.out.bam {output}
        mkdir -p {params.starlogs}
        mv {params.prefix}Log.final.out {params.starlogs}
        mv {params.prefix}Log.out {params.starlogs}
        mv {params.prefix}Log.progress.out {params.starlogs}
        '''
        
rule sort_by_name:
    input: 'mapped/bams/{sample}.bam'
    output: 'mapped/bams/{sample}.sortedByName.bam'
    threads: 4 
    shell:
        '''
        samtools sort -n -T {wildcards.sample} -@ {threads} -o {output} {input}
        '''
        
rule count:
    input: 'mapped/bams/{sample}.sortedByName.bam'
    params:
        phred_cutoff=5
    output: 'mapped/counts_strict/{sample}.counts.tsv'
    shell:
        r'''
        htseq-count --order=name --format=bam --mode=intersection-strict --stranded=no --minaqual={params.phred_cutoff} --type=exon --idattr=gene_id {input} {config[gtf]} > {output}
        '''
        
rule format_counts:
    input: 'mapped/counts_strict/{sample}.counts.tsv'
    output: 'mapped/counts_strict/{sample}.counts.noversion.tsv'
    shell:
        r'''
        cat {input} | sed -E 's/\.[0-9]+//' > {output}
        '''
rule featurecounts:
    input: expand('mapped/bams/{sample}.sortedByName.bam', sample=SAMPLES)
    output: 'mapped/featureCounts/fcounts.tsv'
    threads: 16
    shell:
        r'''featureCounts -a {config[gtf]} -o {output} -t exon -g gene_id -Q 4 -T {threads} {input}'''      
        
rule format_fcounts:
    input: 'mapped/featureCounts/fcounts.tsv'
    output: 'mapped/featureCounts/fcounts.noversion.tsv'
    shell:
        r'''
        cat {input} | sed -E 's/\.[0-9]+//' > {output}

        '''  