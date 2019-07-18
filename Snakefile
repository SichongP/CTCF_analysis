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

print(names)

rule all:
    input: sig = expand("Results/sourmash/{sample}.fq.sig", sample = names),
        bai = expand("Results/mapping/{sample}_sorted.bam.bai", sample = names)

rule sourmash_sig:
    input: "raw_data/{sample}.fq.gz"
    output: "Results/sourmash/{sample}.fq.sig"
    log: "logs/sourmash_sig/{sample}.log"
    benchmark: "benchmarks/sourmash_sig/{sample}.tsv"
    params: time="120"
    conda: "envs/sourmash.yaml"
    shell:
     """
     sourmash --scaled 1000 -K 21,31,51 -o {output} {input}
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
    log: "logs/mapping/bwa_{sample}.log"
    benchmark: "benchmarks/mapping/{sample}.tsv"
    shell:
     """
     bwa mem -t 4 -R "@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tPL:illumina" -o {output} ~/reference/equcab3/equcab3_RefSeq_bwa {input}
     """

rule convert_sam:
    input: "Results/mapping/{sample}.sam"
    output: "Results/mapping/{sample}.bam"
    conda: "envs/samtools.yaml"
    params: time="360"
    log: "logs/mapping/samtools_convert_{sample}.log"
    benchmark: "benchmarks/mapping/convert_{sample}.tsv"
    shell:
     """
     samtools view -hb -@ 2 {input} > {output}
     rm {input}
     """

rule sort_bam:
    input: "Results/mapping/{sample}.bam"
    output: "Results/mapping/{sample}_sorted.bam"
    conda: "envs/samtools.yaml"
    params: time="1200"
    log: "logs/mapping/sort_bam_{sample}.log"
    benchmark: "benchmarks/mapping/sort_bam_{sample}.tsv"
    shell:
     """
     samtools sort -@ 4 -m 2G -o {output} {input}
     rm {input}
     """

rule index_bam:
    input: "Results/mapping/{sample}_sorted.bam"
    output: "Results/mapping/{sample}_sorted.bam.bai"
    conda: "envs/samtools.yaml"
    params: time="600"
    log: "logs/mapping/index_bam_{sample}.log"
    benchmark: "benchmarks/mapping/index_bam_{sample}.tsv"
    shell:
     """
     samtools index -@ 4 {input}
     """
