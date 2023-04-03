# -*- coding: utf-8 -*-
import sys
import argparse


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--project', help='Project identifier', required=True)
    parser.add_argument('-f', '--fasta', help='Contig FASTA file (.fasta/.fasta.gz)', required=True)
    parser.add_argument('-1', '--hic-forward', nargs='+', help='Forward Hi-C reads (.fastq.gz)', required=True)
    parser.add_argument('-2', '--hic-reverse', nargs='+', help='Reverse Hi-C reads (.fastq.gz)', required=True)
    parser.add_argument('-w', '--work-dir', help='Working Directory', required=True)
    parser.add_argument('-g', '--groups', help='Group of Chromosome (Integer)', required=True)
    parser.add_argument('-e', '--enzyme', help='Restriction Enzyme (default: HindIII)', default='HindIII')
    parser.add_argument('-s', '--skip', nargs='*', help='Steps to be skipped')
    return parser.parse_args()


def main():
    args = parse_args()
    # Check whether the project identifier existed or not.
    pass
    return 0


if __name__ == '__main__':
    sys.exit(main())