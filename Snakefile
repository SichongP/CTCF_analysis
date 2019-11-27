localrules: all 
#, collect_stats, install_mspc, getNarrowpeaks, gatherJSD, collect_filtered_stats
ruleorder: filter_bam > convert_sam
import re
import os

configfile: "config.yaml"

# Get file names
## File names should be in the format of: rep_tissue(_Input).fq.gz
files = []
for dirpath, dirname, filenames in os.walk("raw_data/"):
    files.extend(filenames)
    break

names = []
for filename in files:
    name = re.match(r"(.+?)\.(fq|fastq).gz", filename)
    if name:
        names.append(name.group(1))

tissues = config['tissues']
#for file in names:
#    name = re.match(r"(.+?)_Input", file)
#    if name:
#        if name.group(1) not in tissues:
#            tissues.append(name.group(1))

reps = config['reps']

def read_frag_length(wildcards):
    with open(checkpoints.collect_predictd.get(sample=wildcards.sample).output[0]) as file:
        for line in file:
            sample, tag, frag = line.strip().split('\t')
            if sample == wildcards.sample:
                return frag

def get_filtered_read_count(wildcards):
    with open(checkpoints.macs2_get_filterdup_metrics.get(sample=wildcards.sample).output[0]) as file:
        for line in file:
            sample, total, filter = line.strip().split('\t')
            if sample == wildcards.sample:
                return filter

def calc_genome_bg(wildcards):
    read_count = get_filtered_read_count(wildcards)
    frag_len = read_frag_length(wildcards)
    return read_count * frag_len / config['genome_size']

rule all:
    input:
        "Results/metrics/mapping_stats.csv",
        "Results/metrics/filtered_stats.csv",
        "Results/figures/PCA.png",
#        "Results/figures/Fingerprint.png",
        "Results/metrics/QC_metrics.tsv",
        "Results/figures/scatterplot_spearman.png",
        "Results/figures/scatterplot_pearson.png",
        "Results/figures/heatmap_spearman.png",
        "Results/figures/heatmap_pearson.png",
        "Results/metrics/LibQ_Metrics.tsv",
#        expand("Results/macs2/filterdup/{sample}_filterdup.bed", sample = names),
#        expand("Results/macs2/predictd/{sample}_predictedModel.R", sample = tissues),
#        "Results/metrics/predictd_macs2.tsv",
#        "Results/metrics/macs2_filterdup.tsv",
#        expand("Results/macs2/pileup/{sample}_filterdup.pileup.bdg", sample = tissues)
        expand("Results/peaks/{sample}.concensus.bed", sample = tissues),
        expand("Results/mapping/beds/{rep}_{tissue}.bed", rep = reps, tissue = tissues),
        expand("Results/figures/{rep}_{tissue}_Cross_Correlation.pdf", rep = reps, tissue = tissues),
        expand("Results/metrics/FRiP/{rep}_{tissue}_FRiP.txt", rep = reps, tissue = tissues),
        expand("Results/metrics/FRiP/{tissue}_FRiP.txt", tissue = tissues),
        "Results/metrics/peak_summary.tsv",
        expand("Results/peaks_IDR/{tissue}/{tissue}.peaks", tissue = tissues),
        expand("Results/figures/IDR/{tissue}.png", tissue = tissues),
        "Results/metrics/merged_peak_summary.tsv",
#        expand("Results/tagAlign/{rep}_{tissue}{input}.tagAlign", rep = reps, tissue = tissues, input = ["","_Input"]),
        expand("Results/motifs/{tissue}/homerMotifs.all.motifs", tissue = tissues),
        expand("Results/macs2/{rep}_{tissue}_FoldEnrichment.bdg", rep = reps, tissue = tissues),
        expand("Results/IDR/{tissue}_pr1_pooled.tagAlign.gz", rep = reps, tissue = tissues),
        expand("Results/IDR/{tissue}_pr2_pooled.tagAlign.gz", rep = reps, tissue = tissues),
        expand("Results/IDR/{tissue}_Input_pr1_pooled.tagAlign.gz", rep = reps, tissue = tissues),
        expand("Results/IDR/{tissue}_Input_pr2_pooled.tagAlign.gz", rep = reps, tissue = tissues),
        expand("Results/IDR/macs2/{tissue}_{pr}_peaks.narrowPeak", tissue = tissues, pr = ["pr1", "pr2"]),
        expand("Results/IDR/concensus_peaks/mspc/{sample}.concensus.bed", sample = tissues),
        expand("Results/IDR/concensus_peaks/IDR/{sample}/{sample}.peaks", sample = tissues),
        expand("Results/figures/pseudopeaks/idr_{sample}.png", sample = tissues),

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

include: "Rules/Mapping.snakemake"
include: "Rules/MappingStats.snakemake"
include: "Rules/PeakCalling.snakemake"
include: "Rules/QualityMetrics.snakemake"
include: "Rules/PeakMetrics.snakemake"
include: "Rules/Motifs.snakemake"
include: "Rules/IDR.snakemake"
