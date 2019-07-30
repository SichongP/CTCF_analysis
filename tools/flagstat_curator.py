#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 18:27:02 2019

@author: sichong
"""

import os, regex, sys

def read_file(file):
    stats = dict()
    stats['QC_passed'] = dict()
    stats['QC_failed'] = dict()
    try:
        with open(file) as f:
            lines = f.readlines()    
    except Exception as e:
        print(e)
        exit(1)
    m = regex.match("(\d+)\s\+\s(\d+)", lines[0])
    stats['QC_passed']['total'], stats['QC_failed']['total'] = m.captures(1)[0], m.captures(2)[0]
    for i in [1,2,3,5,6,7,9,11,12]:
        m = regex.match("(\d+)\s\+\s(\d+)\s([\w\s\(\)>=\d]+)", lines[i].strip())
        stats['QC_passed'][m.captures(3)[0]], stats['QC_failed'][m.captures(3)[0]] = m.captures(1)[0], m.captures(2)[0]
    # 19846165 + 0 mapped (99.34%:N/A)
    m = regex.match("(\d+)\s\+\s(\d+)\s([\w\s]+)\s\((.+?):(.+?)\)", lines[4])
    stats['QC_passed'][m.captures(3)[0]] = {'count': m.captures(1)[0], 'percent': m.captures(4)[0]}
    stats['QC_failed'][m.captures(3)[0]] = {'count': m.captures(2)[0], 'percent': m.captures(5)[0]}
    m = regex.match("(\d+)\s\+\s(\d+)\s([\w\s]+)\s\((.+?):(.+?)\)", lines[8])
    stats['QC_passed'][m.captures(3)[0]] = {'count': m.captures(1)[0], 'percent': m.captures(4)[0]}
    stats['QC_failed'][m.captures(3)[0]] = {'count': m.captures(2)[0], 'percent': m.captures(5)[0]}
    m = regex.match("(\d+)\s\+\s(\d+)\s([\s\w]+)\s\((.+?):(.+?)\)", lines[10])
    stats['QC_passed'][m.captures(3)[0]] = {'count': m.captures(1)[0], 'percent': m.captures(4)[0]}
    stats['QC_failed'][m.captures(3)[0]] = {'count': m.captures(2)[0], 'percent': m.captures(5)[0]}
    return stats
    
def write_stats(stats, out, name):
    out.write(",".join([name,stats['QC_passed']['total'],stats['QC_failed']['total'],stats['QC_passed']['duplicates'],stats['QC_passed']['mapped']['count'],stats['QC_passed']['mapped']['percent'],stats['QC_passed']['with itself and mate mapped'],stats['QC_passed']['read1'],stats['QC_passed']['read2'],stats['QC_passed']['properly paired']['count'],stats['QC_passed']['properly paired']['percent'],stats['QC_passed']['with itself and mate mapped'],stats['QC_passed']['singletons']['count'],stats['QC_passed']['singletons']['percent'],stats['QC_passed']['with mate mapped to a different chr'],stats['QC_passed']['with mate mapped to a different chr (mapQ>=5)'],stats['QC_failed']['duplicates'],stats['QC_failed']['mapped']['count'],stats['QC_failed']['mapped']['percent'],stats['QC_failed']['with itself and mate mapped'],stats['QC_failed']['read1'],stats['QC_failed']['read2'],stats['QC_failed']['properly paired']['count'],stats['QC_failed']['properly paired']['percent'],stats['QC_failed']['with itself and mate mapped'],stats['QC_failed']['singletons']['count'],stats['QC_failed']['singletons']['percent'],stats['QC_failed']['with mate mapped to a different chr'],stats['QC_failed']['with mate mapped to a different chr (mapQ>=5)'],"\n"]))

def main():
#    in_dir = sys.argv[1]
#    out_path = sys.argv[2]
    in_dir = snakemake.params['dir']
    out_path = snakemake.output[0]
    for dirpath, dirnames, filenames in os.walk(in_dir):
        files = filenames
        path = dirpath
        break
    with open(out_path, 'w') as out:
        out.write("## This is a tabular format curated from sambamba/samtools flagstat output.\n## Author: Sichong Peng (scpeng@ucdavis.edu)\n")
        out.write("File,Total(QC passed),Total(QC failed),Duplicates(QC passed),Mapped(QC passed),Mapped%(QC passed),Paired(QC passed),R1(QC passed),R2(QC passed),Properly paired(QC passed),Properly paired%(QC passed),Paired mapped(QC passed),Mate only(QC passed),Mate only%(QC passed),Mate mapped to diff chr(QC passed), Mate mapped to diff chr(mapQ>5)(QC passed),Duplicates(QC failed),Mapped(QC failed),Mapped%(QC failed),Paired(QC failed),R1(QC failed),R2(QC failed),Properly paired(QC failed),Properly paired%(QC failed),Paired mapped(QC failed),Mate only(QC failed),Mate only%(QC failed),Mate mapped to diff chr(QC failed), Mate mapped to diff chr(mapQ>5)(QC failed)\n")
        for file in files:
            stats = read_file(path+"/"+file)
            write_stats(stats, out, file.rsplit(".",2)[0])

if __name__ == "__main__":
    #pprint.pprint(read_file("data/CP_Adipose.stats.txt"))
    main()
