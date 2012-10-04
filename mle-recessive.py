"""
This is currently a script but later it should be web accessible.
First version is October 4 2012.
It is for testing a pure recessive disease codon preference model
against the genic model of Yang and Nielsen 2008.
This script does not use so many biological libraries
because for some reason I am adopting a more matlab-like style.
"""

import argparse
from itertools import product

import numpy as np
from scipy import optimize, special, linalg

# http://en.wikipedia.org/wiki/Stop_codon
g_stop = {'tag', 'taa', 'tga'}

# http://en.wikipedia.org/wiki/File:Transitions-transversions-v3.png
g_ts = {'ag', 'ga', 'ct', 'tc'}
g_tv = {'ac', 'ca', 'gt', 'tg', 'at', 'ta', 'cg', 'gc'}

# http://en.wikipedia.org/wiki/DNA_codon_table
g_code = {
        ('gct', 'gcc', 'gca', 'gcg'),
        ('cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'),
        ('aat', 'aac'),
        ('gat', 'gac'),
        ('tgt', 'tgc'),
        ('caa', 'cag'),
        ('gaa', 'gag'),
        ('ggt', 'ggc', 'gga', 'ggg'),
        ('cat', 'cac'),
        ('att', 'atc', 'ata'),
        ('tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'),
        ('aaa', 'aag'),
        ('atg',),
        ('ttt', 'ttc'),
        ('cct', 'ccc', 'cca', 'ccg'),
        ('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'),
        ('act', 'acc', 'aca', 'acg'),
        ('tgg',),
        ('tat', 'tac'),
        ('gtt', 'gtc', 'gta', 'gtg'),
        ('taa', 'tga', 'tag'),
        }

def get_fixation_genic(x):
    if not x:
        return 1.0
    else:
        return 1.0 / special.hyp1f1(1.0, 2.0, -x).real

def get_fixation_recessive_disease(x):
    if not x:
        return 1.0
    elif x < 0:
        return 1.0 / special.hyp1f1(0.5, 1.5, -x).real
    elif 0 < x:
        return 1.0 / special.hyp1f1(1.0, 1.5, -x).real

def enum_codons():
    """
    Enumerate lower case codon strings with all three stop codons at the end.
    Speed is not important.
    @return: a list of 64 codons
    """
    codons = [''.join(triple) for triple in product('acgt', repeat=3)]
    return sorted(set(codons) - set(g_stop)) + sorted(g_stop)

def get_ts_tv(codons):
    """
    Get binary matrices defining codon pairs differing by single changes.
    Speed is not important.
    @param codons: sequence of lower case codon strings
    @return: two binary numpy arrays
    """
    ncodons = len(codons)
    ts = np.zeros((ncodons, ncodons), dtype=int)
    tv = np.zeros((ncodons, ncodons), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            nts = sum(1 for p in zip(ci,cj) if ''.join(p) in g_ts)
            ntv = sum(1 for p in zip(ci,cj) if ''.join(p) in g_tv)
            if nts == 1 and ntv == 0:
                ts[i, j] = 1
            if nts == 0 and ntv == 1:
                tv[i, j] = 1
    return ts, tv

def get_syn(codons):
    """
    @return: binary matrix for synonymity
    """


def main(args):
    #
    # get a sequence of codons with stop codons at the end
    codons = enum_codons()
    if len(codons) != 64:
        raise Exception
    #
    # check the genetic code for typos
    table_codons = list(c for cs in g_code for c in cs)
    if len(g_code) != 21:
        raise Exception
    if len(table_codons) != len(set(table_codons)):
        raise Exception
    if set(codons) - set(table_codons):
        raise Exception(set(codons) - set(table_codons))
    if set(table_codons) - set(codons):
        raise Exception(set(table_codons) - set(codons))

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

