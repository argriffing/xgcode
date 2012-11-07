"""
This is for using numpy with codon models.
"""

from itertools import product

import numpy
from numpy import testing


# http://en.wikipedia.org/wiki/Stop_codon
g_stop = {'tag', 'taa', 'tga'}

# http://en.wikipedia.org/wiki/Human_mitochondrial_genetics
g_stop_mito = {'tag', 'taa', 'aga', 'agg'}

# http://en.wikipedia.org/wiki/File:Transitions-transversions-v3.png
g_ts = {'ag', 'ga', 'ct', 'tc'}
g_tv = {'ac', 'ca', 'gt', 'tg', 'at', 'ta', 'cg', 'gc'}

# an ordered tuple of unordered pairs for a gtr model
g_nuc_gtr = ('ac', 'ag', 'at', 'cg', 'ct', 'gt')

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

# http://en.wikipedia.org/wiki/Human_mitochondrial_genetics
g_code_mito = {
        ('gct', 'gcc', 'gca', 'gcg'),
        ('cgt', 'cgc', 'cga', 'cgg'),
        ('aat', 'aac'),
        ('gat', 'gac'),
        ('tgt', 'tgc'),
        ('caa', 'cag'),
        ('gaa', 'gag'),
        ('ggt', 'ggc', 'gga', 'ggg'),
        ('cat', 'cac'),
        ('att', 'atc'),
        ('tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'),
        ('aaa', 'aag'),
        ('atg', 'ata'),
        ('ttt', 'ttc'),
        ('cct', 'ccc', 'cca', 'ccg'),
        ('tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'),
        ('act', 'acc', 'aca', 'acg'),
        ('tgg', 'tga'),
        ('tat', 'tac'),
        ('gtt', 'gtc', 'gta', 'gtg'),
        ('tag', 'taa', 'aga', 'agg'),
        }


def enum_codons(stop):
    """
    Enumerate lower case codon strings with all stop codons at the end.
    @return: a list of 64 codons
    """
    codons = [''.join(triple) for triple in product('acgt', repeat=3)]
    return sorted(set(codons) - set(stop)) + sorted(stop)

def get_hamming(codons):
    """
    Get the hamming distance between codons, in {0, 1, 2, 3}.
    @param codons: sequence of lower case codon strings
    @return: matrix of hamming distances
    """
    ncodons = len(codons)
    ham = numpy.zeros((ncodons, ncodons), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            ham[i, j] = sum(1 for a, b in zip(ci, cj) if a != b)
    return ham

def get_ts_tv(codons):
    """
    Get binary matrices defining codon pairs differing by single changes.
    @param codons: sequence of lower case codon strings
    @return: two binary numpy arrays
    """
    ncodons = len(codons)
    ts = numpy.zeros((ncodons, ncodons), dtype=int)
    tv = numpy.zeros((ncodons, ncodons), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            nts = sum(1 for p in zip(ci,cj) if ''.join(p) in g_ts)
            ntv = sum(1 for p in zip(ci,cj) if ''.join(p) in g_tv)
            if nts == 1 and ntv == 0:
                ts[i, j] = 1
            if nts == 0 and ntv == 1:
                tv[i, j] = 1
    return ts, tv

def get_exch_ts_tv(codons):
    """
    This is a more sophisticated version of get_ts_tv.
    Or alternatively it is a more restricted version of get_gtr.
    It returns an ndim-3 matrix whose shape is (ncodons, ncodons, 2)
    where the third axis specifies transitions and transversions.
    The name exch refers to exchangeability, because this function
    precomputes an ndarray that is used as a component to help build
    the part of the rate matrix that corresponds
    to the nucleotide exchangeability (as opposed to overall rate,
    or nucleotide equilibrium probabilities,
    or mutation-selection codon exchangeability) in the codon rate matrix.
    @param codons: sequence of lower case codon strings
    @return: a numpy array of ndim 3
    """
    ncodons = len(codons)
    ham = get_hamming(codons)
    M = numpy.zeros((ncodons, ncodons, 2), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            if ham[i, j] == 1:
                if any(a+b in g_ts for a,b in zip(ci,cj)):
                    M[i, j, 0] = 1
                if any(a+b in g_tv for a,b in zip(ci,cj)):
                    M[i, j, 1] = 1
    return M

def get_gtr(codons):
    """
    This is a generalization of get_ts_tv_exch.
    It returns a higher dimensional ndarray
    whose shape is (ncodons, ncodons, 6) where the dimension of the last
    axis is the number of upper off-diagonal entries in a nucleotide
    rate matrix, that is, 4*(4-1)/2 = 6.
    The value of M[i, j, k] is 1 if codons i and j differ at exactly
    one nucleotide position and k is the type of the unordered difference,
    otherwise M[i, j, k] is 0.
    This is a very sparse and wasteful representation,
    but it is nice for vectorization.
    @param codons: sequence of lower case codon strings
    @return: a numpy array of ndim 3
    """
    ncodons = len(codons)
    ham = get_hamming(codons)
    M = numpy.zeros((ncodons, ncodons, 6), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            if ham[i, j] == 1:
                for k, pk in enumerate(g_nuc_gtr):
                    if any(ci[a]+cj[a]==pk for a in range(3)):
                        M[i, j, k] = 1
                    if any(cj[a]+ci[a]==pk for a in range(3)):
                        M[i, j, k] = 1
    return M

def get_syn_nonsyn(code, codons):
    """
    Get binary matrices defining synonymous or nonynonymous codon pairs.
    @return: two binary matrices
    """
    ncodons = len(codons)
    inverse_table = dict((c, i) for i, cs in enumerate(code) for c in cs)
    syn = numpy.zeros((ncodons, ncodons), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            if inverse_table[ci] == inverse_table[cj]:
                syn[i, j] = 1
    return syn, 1-syn

def get_compo(codons):
    """
    Get a matrix defining site-independent nucleotide composition of codons.
    @return: integer matrix
    """
    ncodons = len(codons)
    compo = numpy.zeros((ncodons, 4), dtype=int)
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
    @return: a three dimensional matrix
    """
    ncodons = len(codons)
    asym_compo = numpy.zeros((ncodons, ncodons, 4), dtype=int)
    for i, ci in enumerate(codons):
        for j, cj in enumerate(codons):
            for k, nt in enumerate('acgt'):
                asym_compo[i, j, k] = sum(1 for a, b in zip(ci, cj) if (
                    a != b and b == nt))
    return asym_compo

def get_lb_neg_ll(subs_counts):
    """
    Get the lower bound negative log likelihood.
    It uses an entropy-based calculation.
    @param subs_counts: codon substitution counts
    """
    nstates = subs_counts.shape[0]
    counts = []
    for i in range(nstates):
        for j in range(nstates):
            if i < j:
                c = subs_counts[i, j] + subs_counts[j, i]
                counts.append(c)
            elif i == j:
                c = subs_counts[i, j]
                counts.append(c)
    # return the entropy of the unordered pair count vector
    probs = numpy.array(counts, dtype=float) / numpy.sum(counts)
    return -numpy.sum(c*numpy.log(p) for c, p in zip(counts, probs) if c)

def get_lb_expected_subs(ham, subs_counts):
    """
    Get the lower bound of expected substitutions.
    This is expected minimum possible number of substitutions
    between aligned codons.
    """
    return numpy.sum(ham * subs_counts) / float(numpy.sum(subs_counts))



class Test_NumpyCodons(testing.TestCase):

    def _help_test_invariants(self, code, stop):
        all_codons = enum_codons(stop)
        codons = all_codons[:-len(stop)]
        ts, tv = get_ts_tv(codons)
        syn, nonsyn = get_syn_nonsyn(code, codons)
        compo = get_compo(codons)
        asym_compo = get_asym_compo(codons)
        ham = get_hamming(codons)
        gtr = get_gtr(codons)
        # check some invariants
        testing.assert_equal(len(all_codons), 64)
        testing.assert_equal(len(codons), 64 - len(stop))
        testing.assert_equal(numpy.unique(ts), [0, 1])
        testing.assert_equal(numpy.unique(tv), [0, 1])
        testing.assert_equal(numpy.unique(gtr), [0, 1])
        # check the genetic code for typos
        table_codons = list(c for cs in code for c in cs)
        testing.assert_equal(len(code), 21)
        testing.assert_equal(len(table_codons), len(set(table_codons)))
        if set(codons) - set(table_codons):
            raise Exception(set(all_codons) - set(table_codons))
        if set(table_codons) - set(all_codons):
            raise Exception(set(table_codons) - set(all_codons))

    def test_mito_invariants(self):
        self._help_test_invariants(g_code_mito, g_stop_mito)

    def test_plain_invariants(self):
        self._help_test_invariants(g_code, g_stop)

if __name__ == '__main__':
    testing.run_module_suite()

