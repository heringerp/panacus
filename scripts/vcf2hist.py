#!/usr/bin/env python3

import argparse
import logging
import sys
import matplotlib.pyplot as plt
from collections import Counter, defaultdict


def parse_variant(line):
    fields = line.split("\t")
    result = {"chrom": fields[0], "pos": fields[1], "id": fields[2]}
    samples = fields[9:]
    counter = Counter()
    for sample in samples:
        haplotypes = sample.split('|')
        counter.update(haplotypes)
    result["counter"] = counter
    return result


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser(description='Converts a VCF file to a histogram for use in panacus')
    parser.add_argument('filename')

    args = parser.parse_args()
    logging.info(f"Converting {args.filename}")
    no_lines = 0

    cli_call = " ".join(sys.argv)
    allele_counts = defaultdict(int)
    with open(args.filename, 'r') as variants_file:
        for line in variants_file:
            line = line.strip()

            # Skip comments
            if line.startswith('#'):
                continue

            variant = parse_variant(line)
            #print(variant)
            for key, value in variant["counter"].items():
                if key.isdigit():
                    allele = int(key)
                    allele_value = value
                    # Include reference sequence
                    if allele == 0:
                        allele_value += 1
                    allele_counts[allele_value] += 1
            if '0' not in variant["counter"]:
                allele_counts[1] += 1
            no_lines += 1
    allele_counts_list = list(sorted(allele_counts.items()))
    max_count = allele_counts_list[-1][0]

    logging.info(f"Parsed {no_lines} lines")

    print(f"# {cli_call}")
    print("panacus\thist")
    print("count\tallele")
    for i in range(0, max_count + 1):
        print(f"{i}\t{allele_counts[i]}")


if __name__ == '__main__':
    main()
