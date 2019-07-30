#!/usr/bin/env python3
import os, re
def read_file(file):
    with open(file) as f:
        lines = f.readlines()
    for line in lines:
        if "# tag size" in line:
            tag_size = line.split()[11]
        elif "# predicted fragment length" in line:
            frag_len = line.split()[12]
    return tag_size, frag_len

def write_file(samples, out_file):
    with open(out_file, 'w') as o:
        o.write("sample\ttag_size\tfragment_length\n")
        for key in samples:
            o.write(key + "\t" + samples[key][0] + "\t" + samples[key][1] + '\n')
    return True

def main():
    in_dir = snakemake.params['dir']
    out_path = snakemake.output[0]
#    in_dir = "Results/logs/macs2/predictd/"
#    out_path = "test_predictd.tsv"

    for dirpath, dirnames, filenames in os.walk(in_dir):
        files = filenames
        path = dirpath
        break
    samples = {}
    for file in files:
        match = re.match(r"(.+?)_predictd.log", file)
        if match:
            sample = match.group(1)
            tag_size, frag_len = read_file(path+"/"+file)
            samples[sample] = [tag_size, frag_len]
    write_file(samples, out_path)

if __name__ == '__main__':
    main()
