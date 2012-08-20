"""
Reproduce a figure from a publication by Kai Zeng 2010. [UNFINISHED]

The values should be exact.
"""

from StringIO import StringIO
import time

import numpy as np
from scipy import optimize

import Form
import FormOut
import MatrixUtil
import StatsUtil
import kaizeng

def get_form():
    """
    @return: the body of a form
    """
    return []

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    N = 10
    k = 4
    params_list = [
            (0.008, 1, 1, 0, 1.5, 1),
            (0.008, 2, 1, 0, 1.5, 1)]
    allele_histograms = np.zeros((2, N+1))
    for i, params in enumerate(params_list):
        mutation, selection = kaizeng.params_to_mutation_selection(N, params)
        P = kaizeng.get_transition_matrix(N, k, mutation, selection)
        v = MatrixUtil.get_stationary_distribution(P)
        for state_index, counts in enumerate(kaizeng.gen_states(N, k)):
            allele_histograms[i, counts[0]] += v[state_index]
    out = StringIO()
    for hist in allele_histograms:
        h = hist[1:-1]
        h /= np.sum(h)
        for i in range(N-1):
            print >> out, i+1, h[i]
        print >> out
    return out.getvalue()

