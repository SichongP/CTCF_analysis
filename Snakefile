localrules: all, collect_stats
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

tissues = []
for file in names:
    name = re.match(r"(.+?)_Input", file)
    if name:
        if name.group(1) not in tissues:
            tissues.append(name.group(1))

rule all:
    input:
        expand("Results/mapping/{sample}_sorted.bam", sample = names),
        "Results/metrics/mapping_stats.csv",
        expand("Results/macs2/filterdup/{sample}_filterdup.bed", sample = names),
        expand("Results/macs2/predictd/{sample}_predictedModel.R", sample = tissues)

rule sourmash_sig:
    input: "raw_data/{sample}.fq.gz"
    output: "Results/sourmash/{sample}.fq.sig"
    log: "Results/logs/sourmash_sig/{sample}.log"
    benchmark: "Results/benchmarks/sourmash_sig/{sample}.tsv"
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
    log: "Results/logs/sourmash_trim/{sample}.log"
    benchmark: "Results/benchmarks/sourmash_trim/{sample}.tsv"
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
    log: "Results/logs/mapping/bwa_{sample}.log"
    benchmark: "Results/benchmarks/mapping/{sample}.tsv"
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
    log: "Results/logs/mapping/samtools_convert_{sample}.log"
    benchmark: "Results/benchmarks/mapping/convert_{sample}.tsv"
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
    log: "Results/logs/mapping/sort_bam_{sample}.log"
    benchmark: "Results/benchmarks/mapping/sort_bam_{sample}.tsv"
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
    log: "Results/logs/mapping/index_bam_{sample}.log"
    benchmark: "Results/benchmarks/mapping/index_bam_{sample}.tsv"
    shell:
     """
     samtools index -@ {threads} {input}
     """

rule bam_stat:
    input:
        bam = "Results/mapping/{sample}_sorted.bam",
        bai = "Results/mapping/{sample}_sorted.bam.bai"
    output: "Results/mapping/stats/{sample}.stats.txt"
    conda: "envs/sambamba.yaml"
    params: time="60"
    threads: 2
    log: "Results/logs/mapping/stat_bam_{sample}.log"
    benchmark: "Results/benchmarks/mapping/stat_bam_{sample}.tsv"
    shell:
     """
     sambamba flagstat -t {threads} {input} 2>{log} 1>{output}
     """

rule collect_stats:
    input: expand("Results/mapping/stats/{sample}.stats.txt", sample = names)
    params: dir = "Results/mapping/stats/"
    output: "Results/metrics/mapping_stats.csv"
    conda: "envs/stat_curator.yaml"
    script: "tools/flagstat_curator.py"

rule macs2_filterdup:
    input: "Results/mapping/{sample}_sorted.bam"
    output: "Results/macs2/filterdup/{sample}_filterdup.bed"
    params: time="30"
    conda: "envs/macs.yaml"
    log: "Results/logs/macs2/filterdup/{sample}_filterdup.log"
    benchmark: "Results/benchmarks/macs2/{sample}_filterdup.tsv"
    shell:
     """
     macs2 filterdup -g 2497530671 -f BAM -i {input} -o {output} --keep-dup auto 2> {log}
     """

rule macs2_predictd:
    input: "Results/macs2/filterdup/{sample}_filterdup.bed"
    output: "Results/macs2/predictd/{sample}_predictedModel.R"
    log: "Results/logs/macs2/predictd/{sample}_predictd.log"
    conda: "envs/macs.yaml"
    params: time="30"
    shell:
     """
     macs2 predictd -g 2497530671 -i {input} -m 5 50 --rfile {wildcards.sample}_predictedModel.R --outdir Results/macs2/predictd/ 2> {log} 
     """
