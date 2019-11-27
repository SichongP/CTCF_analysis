#!/usr/bin/env python3
# Adopted from Colin Kern
import subprocess

outfile = snakemake.output[0]
tissues = snakemake.params.tissues
reps = snakemake.params.reps
gsize = float(snakemake.config['genome_size'])

with open(outfile, 'w') as out:
    out.write("Tissue\tMethod\tPeaks\tp<0.05\tCoverage\tFRiP_repA\tFRiP_repB\tPR\tPR_p<0.05\tPR_coverage\n")
    for tissue in tissues:
        peakfile = "Results/peaks/{}.concensus.bed_temp".format(tissue)
        output = subprocess.check_output(['wc', '-l', peakfile])
        significant = subprocess.check_output("awk '{{if ($4>1.3) {{print $0}}}}' {} | wc -l".format(peakfile), shell = True)
#        print(peakfile)
        total = sum([abs(int(line.split()[2])-int(line.split()[1])) for line in open(peakfile)])
        out.write('\t'.join([tissue, 'MSPC', str(int(output.split()[0])), str(int(significant.split()[0])), str(total / gsize)]))
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
        peakfile = "Results/IDR/concensus_peaks/mspc/{}.concensus.bed_temp".format(tissue)
        output = subprocess.check_output(['wc', '-l', peakfile])
        significant = subprocess.check_output("awk '{{if ($4>1.3) {{print $0}}}}' {} | wc -l".format(peakfile), shell = True)
#        print(peakfile)
        total = sum([abs(int(line.split()[2])-int(line.split()[1])) for line in open(peakfile)])
        out.write('\t' + '\t'.join([str(int(output.split()[0])), str(int(significant.split()[0])), str(total / gsize)]))
        out.write('\n')
        peakfile = "Results/peaks_IDR/{}/{}.peaks".format(tissue, tissue)
        output = subprocess.check_output(['wc', '-l', peakfile])
        significant = subprocess.check_output("awk '{{if ($5>540) {{print $0}}}}' {} | wc -l".format(peakfile), shell = True)
        total = sum([abs(int(line.split()[2])-int(line.split()[1])) for line in open(peakfile)])
        out.write('\t'.join([tissue, 'IDR', str(int(output.split()[0])), str(int(significant.split()[0])), str(total / gsize)]))
        fripfile = "Results/metrics/FRiP/{}_IDR_FRiP.txt".format(tissue)
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
        peakfile = "Results/IDR/concensus_peaks/IDR/{}/{}.peaks".format(tissue, tissue)
        output = subprocess.check_output(['wc', '-l', peakfile])
        significant = subprocess.check_output("awk '{{if ($5>540) {{print $0}}}}' {} | wc -l".format(peakfile), shell = True)
        total = sum([abs(int(line.split()[2])-int(line.split()[1])) for line in open(peakfile)]) 
        out.write('\t' + '\t'.join([str(int(output.split()[0])), str(int(significant.split()[0])), str(total / gsize)]))
        out.write('\n')

