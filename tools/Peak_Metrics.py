#!/usr/bin/env python3
# Adopted from Colin Kern
import subprocess

outfile = snakemake.output[0]
tissues = snakemake.params.tissues
reps = snakemake.params.reps
gsize = float(snakemake.config['genome_size'])

with open(outfile, 'w') as out:
    out.write("Tissue\tReplicate\tPeaks\tCoverage\tFRiP\n")
    for tissue in tissues:
        for rep in sorted(reps):
            sample = "{}_{}".format(rep, tissue)
            peakfile = "Results/macs2/{}_peaks.narrowPeak".format(sample)
            output = subprocess.check_output(['wc', '-l', peakfile])
            total = sum([abs(int(line.split()[2])-int(line.split()[1])) for line in open(peakfile)])
            with open("Results/metrics/FRiP/{}_FRiP.txt".format(sample)) as f:
                for line in f:
                    cols = line.split()
                    if rep in cols[0]:
                        frip = float(cols[2])
            out.write('\t'.join([tissue, rep, str(int(output.split()[0])), str(total / gsize), str(frip / 100)]) + '\n')
        peakfile = "Results/peaks/{}.concensus.bed_temp".format(tissue)
        output = subprocess.check_output(['wc', '-l', peakfile])
        total = sum([abs(int(line.split()[2])-int(line.split()[1])) for line in open(peakfile)])
        out.write('\t'.join([tissue, 'Combined', str(int(output.split()[0])), str(total / gsize)]))
        fripfile = "Results/metrics/FRiP/{}_FRiP.txt".format(tissue)
        for rep in sorted(reps):
            with open(fripfile) as f:
                for line in f:
                    cols = line.split()
                    if rep in cols[0]:
                        frip = float(cols[2])
                        break
                    else:
                        frip = 0
            out.write('\t' + str(frip / 100))
        out.write('\n')
