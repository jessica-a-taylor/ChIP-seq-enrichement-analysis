#!/usr/bin/env python3

from argparse import ArgumentParser

from Bio import SeqIO

if __name__ == '__main__':
    parser = ArgumentParser(description='Calculate and print total length of FASTA file')

    parser.add_argument('fasta_file', help='Path to the FASTA file')

    args = parser.parse_args()

    total_len = 0
    for seq in SeqIO.parse(args.fasta_file, 'fasta'):
        total_len += len(seq)

    print(total_len)
