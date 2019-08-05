localrules: all, collect_stats, collect_predictd, macs2_get_filterdup_metrics, install_mspc, getNarrowpeaks, gatherJSD
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
        "Results/figures/PCA.png",
        "Results/figures/Fingerprint.png",
        "Results/metrics/QC_metrics.tsv",
#        expand("Results/macs2/filterdup/{sample}_filterdup.bed", sample = names),
#        expand("Results/macs2/predictd/{sample}_predictedModel.R", sample = tissues),
#        "Results/metrics/predictd_macs2.tsv",
#        "Results/metrics/macs2_filterdup.tsv",
#        expand("Results/macs2/pileup/{sample}_filterdup.pileup.bdg", sample = tissues)
        expand("Results/peaks/{sample}.concensus.bed", sample = tissues)
#        expand("Results/metrics/{rep}_{sample}_LibraryQC.tsv", sample = tissues, rep = reps)

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
    output: temp("Results/mapping/{sample}.sam")
    conda: "envs/bwa.yaml"
    params: time="1-0"
    threads: 4
    log: "Results/logs/mapping/bwa_{sample}.log"
    benchmark: "Results/benchmarks/mapping/{sample}.tsv"
    shell:
     """
     bwa mem -t {threads} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:illumina" -o {output} ~/reference/equcab3/equcab3_RefSeq_bwa {input}
     """

rule markDup:
    input: bam = "Results/mapping/{sample}_sorted.bam",
        bai = "Results/mapping/{sample}_sorted.bam.bai"
    output: bam = "Results/mapping/{sample}_markdup.bam",
        bai = "Results/mapping/{sample}_markdup.bam.bai"
    conda: "envs/sambamba.yaml"
    params: time = "120"
    log: "Results/logs/mapping/{sample}_markdup.log"
    benchmark: "Results/benchmarks/mapping/{sample}_markdup.tsv"
    threads: 4
    shell:
     """
     sambamba markdup -t {threads} --tmpdir=/scratch/peng/{wildcards.sample}_markdup/ {input.bam} {output.bam} &> {log}
     rm -rf /scratch/peng/{wildcards.sample}_markdup/
     mv {input.bai} {output.bai}
     """

rule convert_sam:
    input: "Results/mapping/{sample}.sam"
    output: temp("Results/mapping/{sample}.bam")
    conda: "envs/samtools.yaml"
    params: time="360"
    threads: 4
    log: "Results/logs/mapping/samtools_convert_{sample}.log"
    benchmark: "Results/benchmarks/mapping/convert_{sample}.tsv"
    shell:
     """
     samtools view -hb -@ {threads} {input} > {output}
     """

rule sort_bam:
    input: "Results/mapping/{sample}.bam"
    output: temp("Results/mapping/{sample}_sorted.bam")
    conda: "envs/samtools.yaml"
    params: time="180"
    threads: 4
    log: "Results/logs/mapping/sort_bam_{sample}.log"
    benchmark: "Results/benchmarks/mapping/sort_bam_{sample}.tsv"
    shell:
     """
     samtools sort -@ {threads} -m 1G -T /scratch/peng/{sample} -o {output} {input}
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
        bam = "Results/mapping/{sample}_markdup.bam",
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
     macs2 filterdup -g {config[genome_size]} -f BAM -i {input} -o {output} --keep-dup auto 2> {log}
     """

checkpoint macs2_get_filterdup_metrics:
    input: expand("Results/macs2/filterdup/{sample}_filterdup.bed", sample = names)
    output: "Results/metrics/macs2_filterdup.tsv"
    log: "Results/logs/macs2/filterdup/collect.log"
    params: dir = "Results/logs/macs2/filterdup/"
    script: "tools/get_filterdup_metrics.py"

rule macs2_predictd:
    input: "Results/macs2/filterdup/{sample}_filterdup.bed"
    output: "Results/macs2/predictd/{sample}_predictedModel.R"
    log: "Results/logs/macs2/predictd/{sample}_predictd.log"
    conda: "envs/macs.yaml"
    params: time="30"
    shell:
     """
     macs2 predictd -g {config[genome_size]} -i {input} -m 5 50 --rfile {wildcards.sample}_predictedModel.R --outdir Results/macs2/predictd/ 2> {log} 
     """

checkpoint collect_predictd:
    input: expand("Results/macs2/predictd/{rep}_{sample}_predictedModel.R", rep = reps, sample = tissues)
    output: "Results/metrics/predictd_macs2.tsv"
    log: "Results/logs/macs2/predictd/collect.log"
    params: dir = "Results/logs/macs2/predictd/"
    script: "tools/get_tagSize.py"

rule get_pileup_cov:
    input: 
        tsv = "Results/metrics/predictd_macs2.tsv",
        bed = "Results/macs2/filterdup/{sample}_filterdup.bed"
    output: "Results/macs2/pileup/{sample}_filterdup.pileup.bdg"
    log: "Results/logs/macs2/pileup/{sample}_pileup.log"
    benchmark: "Results/benchmarks/macs2/{sample}_pileup.tsv"
    params: time = "30", frag_len = read_frag_length
    conda: "envs/macs.yaml"
    shell:
     """
     macs2 pileup -i {input.bed} -o {output} -f BED --extsize {params.frag_len}
     """

rule slocal_bg:
    input: "Results/macs2/filterdup/{sample}_filterdup.bed"
    output: "Results/macs2/background/{sample}_slocal_bg.bdg"
    log: "Results/logs/macs2/background/{sample}_slocal.log"
    benchmark: "Results/benchmarks/macs2/{sample}_slocal.tsv"
    params: time = "30"
    conda: "envs/macs.yaml"
    shell:
     """
     macs2 pileup -i {input} -B --extsize 500 -o {output} -f BED
     """

rule llocal_bg:
    input: "Results/macs2/filterdup/{sample}_filterdup.bed"
    output: "Results/macs2/background/{sample}_llocal_bg.bdg"
    log: "Results/logs/macs2/background/{sample}_llocal.log"
    benchmark: "Results/benchmarks/macs2/{sample}_llocal.tsv"
    params: time = "30"
    conda: "envs/macs.yaml"
    shell:
     """
     macs2 pileup -i {input} -B --extsize 5000 -o {output} -f BED
     """

rule d_bg:
    input: "Results/macs2/filterdup/{sample}_filterdup.bed"
    output: "Results/macs2/background/{sample}_d_bg.bdg"
    log: "Results/logs/macs2/background/{sample}_d_bg.log"
    benchmark: "Results/benchmarks/macs2/{sample}_d_bg.tsv"
    params: time = "30", frag_len = read_frag_length
    conda: "envs/macs.yaml"
    shell:
     """
     macs2 pileup -i {input} -B -o {output} -f BED --extsize <(echo "{params.frag_len}/2")
     """

rule callpeak:
    input: 
        treatment = "Results/mapping/{sample}_sorted.bam",
        control = "Results/mapping/{sample}_Input_sorted.bam"
    output: "Results/macs2/{sample}_peaks.narrowPeak"
    log: "Results/logs/macs2/{sample}_callpeak.log"
    benchmark: "Results/benchmarks/macs2/{sample}_callpeak.tsv"
    params: time = "60"
    conda: "envs/macs.yaml"
    shell:
     """
     macs2 callpeak -t {input.treatment} -c {input.control} -n {wildcards.sample} --outdir "Results/macs2/" -B --SPMR -g {config[genome_size]} --keep-dup auto 2> {log}
     """

rule getNarrowpeaks:
    input: "Results/macs2/{sample}_peaks.narrowPeak"
    output: temp("Results/peaks/{sample}.bed")
    conda: "envs/blank.yaml"
    shell:
     """
     cut -f 1,2,3,4,8 {input} > {output}
     """

rule install_mspc:
    output: "tools/mspc/dotnet", "tools/mspc/CLI.dll"
    log: "Results/logs/mspc/installation.log"
    conda: "envs/blank.yaml"
    shell:
     """
     bash tools/mspc_install.sh &> {log}
     """

rule concensus_peaks:
    input:
        tool = "tools/mspc/CLI.dll",
        repA = expand("Results/peaks/{rep}_{{sample}}.bed", rep = reps[0]),
        repB = expand("Results/peaks/{rep}_{{sample}}.bed", rep = reps[1])
    output: "Results/peaks/{sample}.concensus.bed"
    log: "Results/logs/mspc/{sample}_concensus_call.log"
    conda: "envs/blank.yaml"
    params: time = "10"
    shell:
     """
     tools/mspc/dotnet tools/mspc/CLI.dll -i {input.repA} -i {input.repB}  -r biological -w 1e-4 -s 1e-8 -a 0.05 -o Results/peaks/{wildcards.sample}/ &> {log}
     cp Results/peaks/{wildcards.sample}/ConsensusPeaks.bed Results/peaks/{wildcards.sample}.concensus.bed
     """

rule multiBigwigSummary:
    input: expand("Results/deeptools/bigwig/{rep}_{sample}_scaled.bw", sample = tissues, rep = reps)
    output: "Results/deeptools/multiBigwigSummary.npz"
    conda: "envs/deeptools.yaml"
    threads: 8
    params: time = "120"
    log: "Results/logs/deeptools/multiBamSummary.log"
    shell:
     """
     multiBigwigSummary bins --bwfiles {input} -o {output} --smartLabels -p {threads} &> {log}
     """

rule bamCompare:
    input: exp = "Results/mapping/{sample}_sorted.bam",
        ctrl = "Results/mapping/{sample}_Input_sorted.bam"
    output: "Results/deeptools/bigwig/{sample}_scaled.bw"
    conda: "envs/deeptools.yaml"
    threads: 6
    params: time = "120"
    log: "Results/logs/deeptools/bamCompare_{sample}.log"
    shell:
     """
     bamCompare -b1 {input.exp} -b2 {input.ctrl} -o {output} -p {threads} --ignoreDuplicates --scaleFactorsMethod SES --effectiveGenomeSize {config[genome_size]} &> {log}
     """

rule plotFingerprint_all:
    input: expand("Results/mapping/{rep}_{sample}_sorted.bam", rep = reps, sample = tissues)
    output: 
        plot = "Results/figures/Fingerprint.png",
    conda: "envs/deeptools.yaml"
    params: time = "360"
    threads: 4
    log: "Results/logs/deeptools/fingerprint.log"
    shell:
     """
     plotFingerprint --bamfiles {input} -o {output.plot} -p {threads} --skipZeros -T "Fingerprint plot" --smartLabels --ignoreDuplicates 
     """

rule calcJSD:
    input: 
        trt = "Results/mapping/{sample}_sorted.bam",
        ctrl = "Results/mapping/{sample}_Input_sorted.bam"
    output: 
        metric = temp("Results/metrics/{sample}_LibraryQC.tsv"),
        plot = "Results/figures/{sample}_Fingerprint.png"
    conda: "envs/deeptools.yaml"
    params: time = "60"
    threads: 4
    log: "Results/logs/deeptools/{sample}_JSD.log"
    shell:
     """
     plotFingerprint --bamfiles {input.trt} {input.ctrl} -o {output.plot} --outQualityMetrics {output.metric} --JSDsample {input.ctrl} -p {threads} --skipZeros -T "Fingerprint plot" --smartLabels --ignoreDuplicates 
     """

rule gatherJSD:
    input: expand("Results/metrics/{rep}_{tissue}_LibraryQC.tsv", rep = reps, tissue = tissues)
    output: "Results/metrics/QC_metrics.tsv"
    script: "tools/gatherJSD.py"
    
rule plotPCA:
    input: "Results/deeptools/multiBigwigSummary.npz"
    output: "Results/figures/PCA.png"
    conda: "envs/deeptools.yaml"
    params: time = "60"
    threads: 4
    log: "Results/logs/deeptools/pca.log"
    shell:
     """
     plotPCA --corData {input} -o {output} -T "PCA plot" --outFileNameData "Results/metrics/pca.tab" --transpose &> {log}
     """
