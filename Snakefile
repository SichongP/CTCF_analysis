import re
import os
# Get file names
files = []
for dirpath, dirname, filenames in os.walk("raw_data/"):
    files.extend(filenames)
    break

names = []
for filename in files:
    name = re.match(r"(.+?)\.(fq|fastq).gz", filename)
    if name:
        names.append(name.group(1))

#print(names)

rule all:
    input: 
        expand("Results/mapping/{sample}_sorted.bam.bai", sample = names),
        expand("Results/mapping/stats/{sample}.stats.txt", sample = names)

rule sourmash_sig:
    input: "raw_data/{sample}.fq.gz"
    output: "Results/sourmash/{sample}.fq.sig"
    log: "logs/sourmash_sig/{sample}.log"
    benchmark: "benchmarks/sourmash_sig/{sample}.tsv"
    params: time="120"
    conda: "envs/sourmash.yaml"
    shell:
     """
     sourmash compute --scaled 1000 -K 21,31,51 -o {output} {input}
     """

rule sourmash_abundance_trim:
    input: "raw_data/{sample}.fq.gz"
    output: "Results/sourmash_trimmed/{sample}.fq.gz.abundtrim"
    conda: "envs/sourmash.yaml"
    params: time="120"
    log: "logs/sourmash_trim/{sample}.log"
    benchmark: "benchmarks/sourmash_trim/{sample}.tsv"
    shell:
     """
     trim-low-abund.py -C 3 -Z 18 -V -M 2e9 -o {output} {input}
     """

rule bwa_mem:
    input: "Results/sourmash_trimmed/{sample}.fq.gz.abundtrim"
    output: "Results/mapping/{sample}.sam"
    conda: "envs/bwa.yaml"
    params: time="1-0"
    threads: 4
    log: "logs/mapping/bwa_{sample}.log"
    benchmark: "benchmarks/mapping/{sample}.tsv"
    shell:
     """
     bwa mem -t {threads} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina" -o {output} ~/reference/equcab3/equcab3_RefSeq_bwa {input}
     """

rule convert_sam:
    input: "Results/mapping/{sample}.sam"
    output: "Results/mapping/{sample}.bam"
    conda: "envs/samtools.yaml"
    params: time="360"
    threads: 4
    log: "logs/mapping/samtools_convert_{sample}.log"
    benchmark: "benchmarks/mapping/convert_{sample}.tsv"
    shell:
     """
     samtools view -hb -@ {threads} {input} > {output}
     rm {input}
     """

rule sort_bam:
    input: "Results/mapping/{sample}.bam"
    output: "Results/mapping/{sample}_sorted.bam"
    conda: "envs/samtools.yaml"
    params: time="180"
    threads: 4
    log: "logs/mapping/sort_bam_{sample}.log"
    benchmark: "benchmarks/mapping/sort_bam_{sample}.tsv"
    shell:
     """
     samtools sort -@ {threads} -m 1G -o {output} {input}
     rm {input}
     """

rule index_bam:
    input: "Results/mapping/{sample}_sorted.bam"
    output: "Results/mapping/{sample}_sorted.bam.bai"
    conda: "envs/samtools.yaml"
    params: time="180"
    threads: 4
    log: "logs/mapping/index_bam_{sample}.log"
    benchmark: "benchmarks/mapping/index_bam_{sample}.tsv"
    shell:
     """
     samtools index -@ {threads} {input}
     """

rule bam_stat:
    input: "Results/mapping/{sample}_sorted.bam"
    output: "Results/mapping/stats/{sample}.stats.txt"
    conda: "envs/sambamba.yaml"
    params: time="60"
    threads: 2
    log: "logs/mapping/stat_bam_{sample}.log"
    benchmark: "benchmarks/mapping/stat_bam_{sample}.tsv"
    shell:
     """
     sambamba flagstat -t {threads} {input} 2>{log} 1>{output}
     """
