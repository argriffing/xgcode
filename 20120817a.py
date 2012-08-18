"""
Reproduce a table from a publication by Kai Zeng 2010. [UNFINISHED]

The table involves sampling so it will not be reproduced exactly.
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

def params_to_mutation_selection(N, params):
    # define the hardcoded number of alleles
    k = 4
    # unpack the params
    theta, ka, kb, g0, g1, g2 = params
    # Expand the parameters into a higher dimensional
    # representation of mutation and selection.
    mutation = np.zeros((k, k))
    for i in range(k):
        for j in range(i+1,k):
            mutation[i,j] = theta / float(2*N)
    for i, j in ((0,1), (0,3), (1,3)):
        mutation[j,i] = ka * mutation[i,j]
    for i, j in ((0,2), (1,2), (2,3)):
        mutation[j,i] = kb * mutation[i,j]
    mutation += np.eye(k) - np.diag(np.sum(mutation, axis=1))
    selection = -np.array([g0, g1, g2, 0]) / float(N)
    return mutation, selection

class G:
    def __init__(self, N, observed_counts):
        """
        @param N: haploid population size
        """
        self.N = N
        self.observed_counts = observed_counts
    def __call__(self, X):
        """
        @param X: six params defining mutation and selection
        @return: negative log likelihood
        """
        # define the hardcoded number of alleles
        k = 4
        # unpack the params
        params = X.tolist()
        theta, ka, kb, g0, g1, g2 = params
        if any(x < 0 for x in (theta, ka, kb)):
            return float('inf')
        mutation, selection = params_to_mutation_selection(self.N, params)
        # get the transition matrix
        P = kaizeng.get_transition_matrix(self.N, k, mutation, selection)
        v = MatrixUtil.get_stationary_distribution(P)
        return -StatsUtil.multinomial_log_pmf(v, self.observed_counts)

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    nsamples = 1
    arr = []
    #
    nsites = 50000
    N = 15*2
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

