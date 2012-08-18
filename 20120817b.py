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
    out = StringIO()
    N = 10
    k = 4
    params_list = [
            (0.008, 1, 1, 0, 1.5, 1),
            (0.008, 2, 1, 0, 1.5, 1)]
    for params in params_list:
        mutation, selection = params_to_mutation_selection(N, params)
        P = kaizeng.get_transition_matrix(N, k, mutation, selection)
        v = MatrixUtil.get_stationary_distribution(P)
    # FIXME
    k = 4
    params = (0.002, 1, 1, 0, 0, 0)
    #params = (0.008, 1, 1, 0.5, 1, 1.5)
    mutation, selection = params_to_mutation_selection(N, params)
    #
    tm = time.time()
    P = kaizeng.get_transition_matrix(N, k, mutation, selection)
    print 'time to construct transition matrix:', time.time() - tm
    #
    tm = time.time()
    v = MatrixUtil.get_stationary_distribution(P)
    print 'time to get stationary distribution:', time.time() - tm
    #
    tm = time.time()
    counts = np.random.multinomial(nsites, v)
    print 'time to sample multinomial counts:', time.time() - tm
    #
    tm = time.time()
    logp = StatsUtil.multinomial_log_pmf(v, counts)
    print 'time to get multinomial log pmf:', time.time() - tm
    #
    for i in range(nsamples):
        counts = np.random.multinomial(nsites, v)
        X0 = np.array(params)
        g = G(N, counts)
        Xopt = optimize.fmin(g, X0)
        arr.append(Xopt)
    print >> out, np.array(arr)
    return out.getvalue()

if __name__ == '__main__':
    k = 4
    nsamples = 100
    settings_list = [
            [15, 50000, [0.002, 1, 1, 0, 0, 0]],
            [15, 50000, [0.002, 1, 2, 0.4, -1.2, 2]],
            [10, 10000, [0.008, 1, 1, 0.5, 1, 1.5]],
            [6, 5000, [0.01, 1, 2, 0, 0, 0]]]
    for N_diploid, nsites, params in settings_list:
        N = 2*N_diploid
        print 'diploid population size = %s, sequence length = %s' % (
                N_diploid, nsites)
        print '\t'.join(str(x) for x in ['Input'] + params)
        mutation, selection = params_to_mutation_selection(N, params)
        P = kaizeng.get_transition_matrix(N, k, mutation, selection)
        v = MatrixUtil.get_stationary_distribution(P)
        arr = []
        for i in range(nsamples):
            counts = np.random.multinomial(nsites, v)
            X0 = np.array(params)
            g = G(N, counts)
            Xopt = optimize.fmin(g, X0)
            arr.append(Xopt.tolist())
        means = []
        cis = []
        for mles in zip(*arr):
            means.append(np.mean(mles))
            x = sorted(mles)
            cis.append([x[2], x[-2]])
        print '\t'.join(str(x) for x in ['Mean MLE'] + means)
        print '\t'.join(str(x) for x in ['95\% CI'] + cis)
        print

