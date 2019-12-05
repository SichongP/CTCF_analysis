#!/usr/bin/env python3

import re, os

files = snakemake.input
out = snakemake.output[0]

with open(out, 'w') as o:
    o.write("Tissue\tComparison\tintersection\tunion\tjaccard\tn_intersections\n")
    for file in files:
        with open(file) as f:
            lines = f.readlines()
        intersection, union, jaccard, n_intersections = lines[1].split()
        filename = os.path.basename(file)
        match = re.match('(.+?)_(.+?)\.txt', filename)
        comparison = match[2]
        tissue = match[1]
        o.write('\t'.join([tissue, comparison, str(intersection), str(union), str(jaccard), str(n_intersections)]) + '\n')

