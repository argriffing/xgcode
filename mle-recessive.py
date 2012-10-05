"""
This is currently a script but later it should be web accessible.
First version is October 4 2012.
It is for testing a pure recessive disease codon preference model
against the genic model of Yang and Nielsen 2008.
This script does not use so many biological libraries
because for some reason I am adopting a more matlab-like style.
All matrices in this script are ndarrays as opposed to actual python matrices.
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
    """
    This is an h function in the notation of Yang and Nielsen.
    This should be applicable entrywise to an ndarray.
    Speed is important.
    """
    return 1.0 / special.hyp1f1(1.0, 2.0, -x).real

def get_fixation_recessive_disease(x):
    """
    This is an h function in the notation of Yang and Nielsen.
    This should be applicable entrywise to an ndarray.
    Speed is important.
    """
    xneg = np.clip(x, -np.inf, 0)
    xpos = np.clip(x, 0, np.inf)
    a = 1.0 / special.hyp1f1(0.5, 1.5, -xneg).real
    b = 1.0 / special.hyp1f1(1.0, 1.5, -xpos).real
    return a + b - 1.0

def enum_codons():
    """
    Enumerate lower case codon strings with all three stop codons at the end.
    Speed does not matter.
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

def get_syn_nonsyn(codons):
    """
    Get binary matrices defining synonymous or nonynonymous codon pairs.
    Speed is not important.
    @return: two binary matrices
    """
    ncodons = len(codons)
    inverse_table = dict((c, i) for i, cs in enumerate(g_code) for c in cs)
    syn = np.zeros((ncodons, ncodons), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            if inverse_table[ci] == inverse_table[cj]:
                syn[i, j] = 1
    return syn, 1-syn

def get_compo(codons):
    """
    Get a matrix defining site-independent nucleotide composition of codons.
    Speed is not important.
    @return: integer matrix
    """
    compo = np.zeros((ncodons, 4), dtype=int)
    for i, c in enumerate(codons):
        for j, nt in enumerate('acgt'):
            compo[i, j] = c.count(nt)
    return compo

def get_asym_compo(codons):
    """
    This is an ugly function.
    Its purpose is to help determine which is the mutant nucleotide type
    given an ordered pair of background and mutant codons.
    This determination is necessary if we want to follow
    the mutation model of Yang and Nielsen 2008.
    Entry [i, j, k] of the returned matrix gives the number of positions
    for which the nucleotides are different between codons i and j and
    the nucleotide type of codon j is 'acgt'[k].
    Speed is not important.
    @return: a three dimensional matrix
    """
    asym_compo = np.zeros((ncodons, ncodons, 4), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            for k, nt in enumerate('acgt'):
                asym_compo[i, j, k] = sum(1 for a, b in zip(ci, cj) if (
                    a != b and b == nt))
    return asym_compo

def get_selection_F(log_counts, compo, log_nt_weights):
    """
    The F and S notation is from Yang and Nielsen 2008.
    Note that three of the four log nt weights are free parameters.
    One of the four log weights is zero and the other three
    are free parameters to be estimated jointly in the
    maximimum likelihood search,
    so this function is inside the optimization loop.
    Speed matters.
    @param log_counts: logs of empirical codon counts
    @param compo: codon composition as defined in the get_compo function
    @param log_nt_weights: un-normalized log mutation process probabilities
    @return: a log selection for each codon, up to an additive constant
    """
    return log_counts - np.dot(compo, log_nt_weights)

def get_selection_S(F):
    """
    The F and S notation is from Yang and Nielsen 2008.
    Speed matters.
    @param F: a selection value for each codon, up to an additive constant
    @return: selection differences F_j - F_i, also known as S_ij
    """
    e = np.ones_like(F)
    return np.outer(e, F) - np.outer(F, e)

def get_fixation(h, S):
    """
    Notation is from Yang and Nielsen 2008.
    Speed matters.
    @param h: apply this to each element of S
    @param S: selection S_ij
    @param: reversible fixation rates with log(F) stationary process
    """
    np.zeros
    pass

def get_Q(
        ts, tv, syn, nonsyn, asym_compo,
        log_mu, log_kappa, log_omega, log_nt_weights, fixation):
    """
    Notation is from Yang and Nielsen 2008.
    Speed matters.
    @param ts: indicator for transition
    @param tv: indicator for transversion
    @param syn: indicator for synonymous codons
    @param nonsyn: indicator for nonsynonymous codons
    @param asym_compo: tensor from get_asym_compo function
    @param log_mu: free param for scaling
    @param log_kappa: free param for transition transversion rate distinction
    @param log_omega: free param for syn nonsyn rate distinction
    @param log_nt_weights: mostly free param array for mutation equilibrium
    @param fixation: matrix of h(S_ij) = h(F_j - F_i)
    @return: rate matrix Q
    """
    mu = math.exp(log_mu)
    kappa = math.exp(log_kappa)
    omega = math.exp(log_omega)
    return mu * (kappa * ts + tv) * (omega * nonsyn + syn) * np.exp(np.dot(
            asym_compo, log_nt_weights)) * fixation

def main(args):
    #
    # precompute some ndarrays
    codons = enum_codons()
    #
    # check invariants of precomputed ndarrays
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
    #
    # Check reversibility of h functions with respect to F,
    # in the notation of Yang and Nielsen 2008.
    for h in (get_fixation_genic, get_fixation_recessive_disease):
        F = np.array([1.2, 2.3, 0, -1.1])
        S = get_selection_S(F)
        fixation = h(S)
        log_ratio = np.log(fixation / fixation.T)
        if not np.allclose(S, log_ratio):
            raise Exception((S, log_ratio))

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

