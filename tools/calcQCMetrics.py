#!/usr/bin/env python3
import subprocess

outfile = snakemake.output[0]

samples = snakemake.params.samples

with open(outfile, 'w') as out:
    out.write("Tissue\tRep\tDepth\tNRF\tPBC1\tPBC2\tNSC\tRSC\tJSD\tSJSD\tIJSD\n")
    for sample in samples:
        rep = sample.split('_')[0]
        tissue = sample.split('_')[1]
        output = [tissue, rep]
        depth = int(subprocess.check_output(["samtools","view","-c","Results/mapping/{}_filtered.bam".format(sample)]))
        result = subprocess.check_output(["bamToBed -i Results/mapping/{}_markdup.bam | awk 'BEGIN{{OFS=\"\t\"}}{{print $1,$2,$3,$6}}' | sort | uniq -c | awk 'BEGIN{{mt=1;m0=1;m1=1;m2=1}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}}'".format(sample)], shell = True)
        metrics = result.split()
        output.extend([str(depth),str(float(metrics[4])),str(float(metrics[5])),str(float(metrics[6]))])
        with open("Results/metrics/spp/{}_spp_stats.txt".format(sample)) as sppf:
            sppstats = sppf.read().split()
        output.extend([str(sppstats[8]), str(sppstats[9])])
        with open("Results/metrics/{}_LibraryQC.tsv".format(sample)) as f:
            lines = f.readlines()
        chipstats = lines[1].split()
        try:
            sjsd = float(chipstats[8])
        except IndexError:
            sjsd = 0
        if len(lines) > 2:
            jsd = float(chipstats[7])
            ijsd = float(lines[2].split()[8])
            if ijsd > sjsd:
                jsd = -jsd
            output.extend([str(jsd), str(sjsd), str(ijsd)])
        else:
            output.extend(["N/A",sjsd,"N/A"])
        out.write('\t'.join(output))
        out.write('\n')
