"""
Extract a 4x4 nucleotide subsitution count matrix.

Get it from the file availble from Ziheng Yang.
It is the mitochondrial data file from Yang-Nielsen 2008.
IMPORTANT: note that the ordering is AGCT not ACGT
"""

import argparse

import numpy as np

g_mtdna_names = {'human_horai', 'chimp_horai'}


def read_yang_mtdna_alignment(lines):
    name = None
    alignment = []
    segments = []
    for line in lines:
        line = line.strip().lower()
        if line in g_mtdna_names or not line:
            if segments:
                dna = ''.join(segments)
                alignment.append((name, dna))
        if line in g_mtdna_names:
            name = line
            segments = []
        elif not line:
            name = None
        else:
            if name:
                segments.append(line)
    return alignment

def main(args):
    with open(args.infile) as fin:
        alignment = read_yang_mtdna_alignment(fin)
    nt_to_i = dict((nt, i) for i, nt in enumerate('agct'))
    M = np.zeros((4, 4), dtype=int)
    d = dict(alignment)
    human_mtdna = d['human_horai']
    chimp_mtdna = d['chimp_horai']
    for a, b in zip(human_mtdna, chimp_mtdna):
        M[nt_to_i[a], nt_to_i[b]] += 1
    print M

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--infile',
            help='codon alignment input file using the format of Ziheng Yang',
            required=True)
    main(parser.parse_args())

