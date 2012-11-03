"""
This script transforms the format of a Yang-Nielsen 2008 data set.

The data set is non-mtDNA and I want to transform it into the format
used by the ntDNA data set from that paper.
The data is available at
http://abacus.gene.ucl.ac.uk/ziheng/data.html
on the website of Ziheng Yang.
The idea is that there are instructions in PAML on how to reproduce
table 1 in that paper, but I want to reproduce table 2 which has
a different format.
My approach is to transform the table 2 formatted data into the format
of the table 1 data and then use the PAML instructions for table 1.
"""

import math
import argparse

import npcodon

def read_yang_alignments(lines):
    """
    Read the file that Ziheng Yang put on his website.
    @param lines: text lines of the alignment
    @return: list of alignments
    """
    names = {'hg18', 'pantro2', 'rhemac2', 'mm8', 'rn4'}
    # read the alignments
    alignments = []
    alignment = []
    expected_len = None
    for line in lines:
        line = line.strip().lower()
        if not line:
            continue
        elements = line.split()
        if elements[0] in names:
            seq = elements[1:]
            if len(seq) * 3 != expected_len:
                raise Exception((len(seq) * 3, expected_len))
            alignment.append((elements[0], seq))
        elif int(elements[0]) == len(names):
            if len(elements) != 2:
                raise Exception(elements)
            expected_len = int(elements[1])
            if alignment:
                alignments.append(alignment)
            alignment = []
        else:
            raise Exception(elements[0])
    if alignment:
        alignments.append(alignment)
    return alignments


def main(args):
    t1, t2 = args.t1, args.t2
    with open(args.infile) as fin:
        alignments = read_yang_alignments(fin)
    arr1 = []
    arr2 = []
    for alignment in alignments:
        d = dict(alignment)
        arr1.append(''.join(d[t1]))
        arr2.append(''.join(d[t2]))
    s1 = ''.join(arr1).upper()
    s2 = ''.join(arr2).upper()
    n = len(s1)
    #print len(s1)
    #print len(s2)
    #
    # this is the transformed file
    print '  2', n
    print
    print t1
    for i in range(int(math.ceil(n/60.0))):
        print s1[60*i : 60*(i+1)]
    print
    print t2
    for i in range(int(math.ceil(n/60.0))):
        print s2[60*i : 60*(i+1)]
    print

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--infile',
            help='codon alignment input file using the format of Ziheng Yang',
            required=True)
    parser.add_argument(
            '--t1',
            default='mm8',
            help='name of first taxon')
    parser.add_argument(
            '--t2',
            default='rn4',
            help='name of second taxon')
    main(parser.parse_args())
