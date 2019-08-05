#!/usr/bin/env python3

files = snakemake.input
out = snakemake.output[0]

def read_file(file):
    with open(file) as f:
        lines = f.readlines()
    return lines

def main():
    with open(out, 'w') as o:
        o.write(read_file(files[0])[0])
        for file in files:
            o.write(read_file(file)[1])

if __name__ == '__main__':
    main()
