#!/usr/bin/env python3

import argparse
import logging
import sys
from collections import Counter, defaultdict


def parse_variant(line, filter):
    fields = line.split("\t")
    result = {"chrom": fields[0], "pos": fields[1], "id": fields[2]}
    samples = fields[9:]
    counter = Counter()
    no_haplotypes = 0
    if filter is not None:
        filter_fields = filter.split("=")
        filter_key = filter_fields[0]
        filter_value = filter_fields[1]
        infos = fields[7].split(";")
        info_values = [info.split("=") for info in infos]
        filter_info = [info[1] for info in info_values if info[0] == filter_key]
        if len(filter_info) >= 1:
            filter_info = filter_info[0]
            info_alleles = [el == filter_value for el in filter_info.split(",")]
            valid_alleles = [str(el + 1) for el in range(0, len(info_alleles)) if info_alleles[el]]
        else:
            valid_alleles = []
    for sample in samples:
        haplotypes = sample.split('|')
        no_haplotypes += len(haplotypes)
        if filter is not None:
            haplotypes = [h if h in valid_alleles else "0" for h in haplotypes]
        counter.update(haplotypes)
    del counter['0']
    result["counter"] = counter
    return result, no_haplotypes


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser(description='Converts a VCF file to a histogram for use in panacus')
    parser.add_argument('filename')
    parser.add_argument('-f', '--filter')

    args = parser.parse_args()
    logging.info(f"Converting {args.filename}")
    no_lines = 0
    no_haplotypes = -1
    filter = args.filter

    cli_call = " ".join(sys.argv)
    allele_counts = defaultdict(int)
    with open(args.filename, 'r') as variants_file:
        for line in variants_file:
            line = line.strip()

            # Skip comments
            if line.startswith('#'):
                continue

            variant, curr_no_haplotypes = parse_variant(line, filter)
            if no_haplotypes == -1:
                no_haplotypes = curr_no_haplotypes + 1
            # print(variant)
            for key, value in variant["counter"].items():
                if key.isdigit():
                    allele = int(key)
                    allele_value = value
                    # Include reference sequence
                    if allele == 0:
                        allele_value += 1
                    allele_counts[allele_value] += 1
            if '0' not in variant["counter"] and False:
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
    for i in range(max_count + 1, no_haplotypes + 1):
        print(f"{i}\t0")


if __name__ == '__main__':
    main()
