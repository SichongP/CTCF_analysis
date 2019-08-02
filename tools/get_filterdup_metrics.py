#!/usr/bin/env python3
import os, re
def read_file(file):
    with open(file) as f:
        lines = f.readlines()
    for line in lines:
        if "tags after filtering" in line:
            filtered_count = line.split()[13]
        elif "total tags in alignment" in line:
            total_count = line.split()[12]
    return total_count, filtered_count

def write_file(samples, out_file):
    with open(out_file, 'w') as o:
        o.write("sample\tTotal_read_count\tAfter_filter\n")
        for key in samples:
            o.write(key + "\t" + samples[key][0] + "\t" + samples[key][1] + '\n')
    return True

def main():
    in_dir = snakemake.params['dir']
    out_path = snakemake.output[0]
#    in_dir = "Results/logs/macs2/filterdup/"
#    out_path = "test_filterdup.tsv"

    for dirpath, dirnames, filenames in os.walk(in_dir):
        files = filenames
        path = dirpath
        break
    samples = {}
    for file in files:
        match = re.match(r"(.+?)_filterdup.log", file)
        if match:
            sample = match.group(1)
            total, filter = read_file(path+"/"+file)
            samples[sample] = [total, filter]
    write_file(samples, out_path)

if __name__ == '__main__':
    main()
