"""
Check a spectral tree reconstruction that uses SVD for branch lengths.

Check that the Laplacian properties are preserved over the iterations.
"""

from StringIO import StringIO
import random
import itertools

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import combobreaker
import nhj

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ntaxa', 'vertices in the completely connected graph',
                5, low=2, high=10),
            Form.Integer('nsamples', 'max iterations', 10000, low=1),
            Form.RadioGroup('split_criterion', 'split criterion', [
                Form.RadioItem('split_fiedler', 'Fiedler', True),
                Form.RadioItem('split_random', 'random')]),
            ]
    return form_objects

def get_form_out():
    return FormOut.Report('summary')

def sample_laplacian_matrix(n):
    """
    @param n: number of states
    @return: a random n x n weighted combinatorial laplacian matrix
    """
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            x = random.expovariate(1)
            L[i, j] = -x
            L[j, i] = -x
    L -= np.diag(np.sum(L, axis=1))
    return L

def foo():
    pass

class Accumulator:
    def __init__(self, ntaxa, fsplit):
        self.ntaxa = ntaxa
        self.fsplit = fsplit
        self.rows = ['hello world']
        self.counterexample = 'no counterexample'
    def __call__(self):
        n = self.ntaxa
        # sample a laplacian matrix
        L = sample_laplacian_matrix(n)
        # get the new vertex and the two new replacement supervertices
        v_gen = itertools.count(n)
        sv_gen = itertools.count(1)
        v_new = next(v_gen)
        sv_new_a = next(sv_gen)
        sv_new_b = next(sv_gen)
        # define the initial structures
        v_to_svs = dict((v, set([0])) for v in range(n))
        sv_to_vs = {0 : set(range(n))}
        edge_to_weight = {}
        for pair in itertools.combinations(range(n), 2):
            edge_to_weight[frozenset(pair)] = random.expovariate(1)
        sv_heap = [(-n, 0)]
        # do the split
        nhj.split(
                v_to_svs, sv_to_vs, edge_to_weight, sv_heap,
                v_new, sv_new_a, sv_new_b, self.fsplit)
        # accumulate
        #self.rows.append('hello world')
    def __str__(self):
        out = StringIO()
        print >> out, '\n'.join(str(row) for row in self.rows)
        print >> out, self.counterexample
        return out.getvalue()

def get_response_content(fs):
    niterations = fs.nsamples
    ntaxa = fs.ntaxa
    nseconds = 4
    if fs.split_fiedler:
        fsplit = nhj.fsplit_fiedler
    elif fs.split_random:
        fsplit = nhj.fsplit_random
    accum = Accumulator(ntaxa, fsplit)
    info = combobreaker.run_callable(
            accum, nseconds=nseconds, niterations=niterations)
    out = StringIO()
    print >> out, str(info)
    return out.getvalue()

