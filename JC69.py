"""
Functions dealing with the JC69 nucleotide rate matrix.

This is the continuous time Markov model defined by Jukes and Cantor in 1969.
Some analysis was attempted using Theano but it is very slow
and is probably not useful except for testing.
"""

from StringIO import StringIO
import unittest
import random
import math
from math import log, exp, expm1

import numpy as np
import theano
import theano.tensor as T

import MatrixUtil
import Util

class FisherInformationTheano:
    def __init__(self, mu, N=4):
        """
        This is the Fisher information for time.
        Try using Theano with the naive definition of Fisher information.
        @param mu: randomization rate
        @param N: number of states
        """
        self.mu = mu
        self.N = N
        self.cached_theano_f = None
        self.cached_theano_df = None
    def _get_theano_f_df(self):
        if not self.cached_theano_f:
            f, df = self._compute_theano_f_df()
            self.cached_theano_f = f
            self.cached_theano_df = df
        return self.cached_theano_f, self.cached_theano_df
    def _compute_theano_f_df(self):
        mu = T.dscalar('mu')
        N = T.dscalar('N')
        t = T.dscalar('t')
        x = T.exp(-t*mu)
        Pii = x + (1-x)/N
        Pij = (1-x)/N
        Jii = Pii/N
        Jij = Pij/N
        ii_contrib = Jii*T.grad(T.log(Jii), t)**2
        ij_contrib = Jij*T.grad(T.log(Jij), t)**2
        fi = N*(N-1)*ij_contrib + N*ii_contrib
        fi_grad = T.grad(fi, t)
        #print theano.pp(fi)
        #print theano.pp(fi_grad)
        f = theano.function([mu, N, t], fi)
        df = theano.function([mu, N, t], fi_grad)
        #print theano.pp(f.maker.env.outputs[0])
        #print theano.pp(df.maker.env.outputs[0])
        return f, df
    def __call__(self, t):
        f, df = self._get_theano_f_df()
        return f(self.mu, self.N, t)
    def deriv(self, t):
        f, df = self._get_theano_f_df()
        return df(self.mu, self.N, t)

class FisherInformation:
    def __init__(self, mu, N=4):
        """
        This is the Fisher information for time.
        Try using Theano with the naive definition of Fisher information.
        @param mu: randomization rate
        @param N: number of states
        """
        self.mu = mu
        self.N = N
    def deriv(self, t):
        """
        This is from a symbolic differentiation output.
        The following expression was differentiated with respect to x.
        N*(N-1)*(1/N)*(N/(1-x))*(u*x/N)**2 +
        N*(1/N)*(1/(x + (1-x)/N))*(-u*x + (u*x)/N)**2
        Alternate form
        - ( (N-1)*(mu*x)**2 ) / ( (x-1)*((N-1)*x+1) )
        """
        x = exp(-self.mu*t)
        N = self.N
        mu = self.mu
        a = (N-1)*mu*mu*x*((N-2)*x+2)
        b = ( (x-1)*((N-1)*x+1) )**2
        return -mu*x*(a/b)
    def __call__(self, t):
        x = exp(-self.mu*t)
        N = self.N
        mu = self.mu
        Pii = x + (1-x)/N
        Pij = (1-x)/N
        Pii_grad = -mu*x + (mu/N)*x
        Pij_grad = (mu/N)*x
        ii_contrib = (Pii_grad**2) / (N*Pii)
        ij_contrib = (Pij_grad**2) / (N*Pij)
        return N*(N-1)*ij_contrib + N*ii_contrib



class MutualInformation:
    def __init__(self, mu, N=4):
        """
        This is the mutual information between X(0) and X(t).
        @param mu: randomization rate
        @param N: number of states
        """
        self.mu = mu
        self.N = N
        # define the logical entropy of the stationary distribution
        self.h = (N-1) / float(N)
    def deriv(self, t):
        x = exp(-self.mu*t)
        return self.h*(log(1-x) - log(x*(self.N-1) + 1))*self.mu*x
    def __call__(self, t):
        x = exp(-self.mu*t)
        xcomp = -expm1(-self.mu*t)
        if not x:
            return 0
        elif not xcomp:
            shannon_entropy = log(self.N)
            return shannon_entropy
        try:
            a = self.h*xcomp*log(xcomp)
            b = (x + xcomp/self.N)*log(self.N*x + xcomp)
            mi = a + b
        except ValueError as e:
            print 'mu:', self.mu
            print 'N:', self.N
            print 'x:', x
            print '1-x:', xcomp
            raise e
        return mi

class IdentitySlopeInformation:
    """
    This is like an information.
    """
    def __init__(self, mu, N=4):
        """
        @param mu: randomization rate
        @param N: number of states
        """
        self.f_identity = ExpectedIdentity(mu, N)
    def deriv(self, t):
        return -self.f_identity.d2(t)
    def __call__(self, t):
        return -self.f_identity.deriv(t)

class ExpectedIdentity:
    """
    The negative derivative of this function is like an information.
    """
    def __init__(self, mu, N=4):
        """
        This is P(X(0) == X(t)) for N-state Jukes-Cantor.
        @param mu: randomization rate
        @param N: number of states
        """
        self.mu = mu
        self.N = N
        # define the logical entropy of the stationary distribution
        self.h = (N-1) / float(N)
    def d2(self, t):
        return self.h * self.mu * self.mu * math.exp(-self.mu * t)
    def deriv(self, t):
        return -self.h * self.mu * math.exp(-self.mu * t)
    def inv(self, p):
        return -math.log((p + self.h - 1) / self.h) / self.mu
    def __call__(self, t):
        return self.h * math.exp(-self.mu * t) + (1 - self.h)

def distance_to_probability(d):
    """
    @param d: the expected number of substitutions along a branch
    @return: the probability that the endpoints are different
    """
    return (3.0 / 4.0) * (1.0 - math.exp(-(4.0 / 3.0) * d))

def probability_to_distance(p):
    """
    Note that p can also be an observed proportion of nucleotide differences.
    In this case the returned value is the maximum likelihood estimate
    of the expected number of substitutions along the branch.
    @param p: the probability that the endpoints of a branch are different
    @return: the expected number of substitutions along the branch
    """
    if p >= (3.0 / 4.0):
        return float('inf')
    return -(3.0 / 4.0) * math.log(1.0 - (4.0 / 3.0) * p)

def sample_xtree_sequences(root, sequence_length):
    """
    @param root: the root of a phylogenetic xtree with branch lengths
    @param sequence_length: the number of nucleotides per sequence
    @return: a list of nucleotide sequences ordered by the xtree leaf labels
    """
    ntips = len(root.get_labeled_vertices())
    sequences = [None]*ntips
    id_to_sequence = {}
    for vertex in root.get_preorder_vertices():
        if vertex.branch:
            letters = []
            p_randomize = 1 - math.exp(-(4.0/3.0)*vertex.branch.length)
            for parent_nt in id_to_sequence[id(vertex.branch.parent)]:
                if random.random() < p_randomize:
                    letters.append(random.choice('ACGT'))
                else:
                    letters.append(parent_nt)
        else:
            letters = [random.choice('ACGT') for i in range(sequence_length)]
        sequence = ''.join(letters)
        if vertex.has_label():
            sequences[vertex.label] = sequence
        else:
            id_to_sequence[id(vertex)] = sequence
    return sequences

def sample_sequences(tree, order, sequence_length):
    """
    This sampler is supposed to be somewhat fast.
    @param tree: a newick tree
    @param order: the requested order of the taxa
    @param sequence_length: the requested length of the sequences
    @return: a list of nucleotide sequences in the requested order
    """
    # assert that the taxa named in the order are the same as the names of the tips
    sorted_tree_tip_names = tuple(sorted(node.get_name() for node in tree.gen_tips()))
    sorted_order_names = tuple(sorted(order))
    assert sorted_tree_tip_names == sorted_order_names
    # get the number of nodes in the tree, including internal nodes
    n = len(list(tree.preorder()))
    # make a preorder array representation of the tree, setting the defaults for the root
    index_to_parent_index = [-1]*n
    index_to_branch_length = [-1]*n
    # map the node id and the node name to the index
    name_to_index = {}
    id_to_index = {}
    for i, node in enumerate(tree.preorder()):
        name_to_index[node.get_name()] = i
        id_to_index[id(node)] = i
        if i:
            index_to_parent_index[i] = id_to_index[id(node.get_parent())]
            index_to_branch_length[i] = node.get_branch_length()
    # for each branch calculate the probability that four way randomization is not needed
    index_to_pequal = [-1]*n
    for i in range(1, n):
        distance = index_to_branch_length[i]
        index_to_pequal[i] = math.exp(-(4.0/3.0)*distance)
    # simulate nucleotide columns
    columns = []
    for iteration in range(sequence_length):
        # choose the state at the root according to the JC69 stationary distribution
        column = [random.choice('ACGT')]
        # choose the rest of the states depending on the state of the parent
        for i in range(1, n):
            parent_state = column[index_to_parent_index[i]]
            if random.random() < index_to_pequal[i]:
                column.append(parent_state)
            else:
                column.append(random.choice('ACGT'))
        columns.append(column)
    # convert the columns to sequences
    sequences = zip(*columns)
    # return the list of sequences in the correct order
    ordered_sequences = [''.join(sequences[name_to_index[name]]) for name in order]
    return ordered_sequences

def get_ML_distance(sa, sb):
    """
    Use a closed form maximum likelihood estimator.
    @param sa: a nucleotide sequence
    @param sb: a nucleotide sequence
    @return: the estimated expected number of changes per position
    """
    assert len(sa) == len(sb)
    assert set(sa+sb) <= set('ACGT')
    n = len(sa)
    n2 = Util.hamming_distance(sa, sb)
    n1 = n - n2
    if n2 >= 3 * n1:
        return float('inf')
    mle = -(3.0/4.0) * math.log((3.0*n1 - n2) / (3.0*n1 + 3.0*n2))
    return mle

def get_ML_distance_matrix(sequence_list):
    """
    Use a closed form maximum likelihood estimator.
    @param sequence_list: an ordered list of nucleotide sequences
    @return: a row major distance matrix
    """
    return MatrixUtil.list_to_matrix(sequence_list, get_ML_distance)


class TestJC69(unittest.TestCase):

    def test_conversion_special_cases(self):
        d = probability_to_distance(0.0)
        self.assertAlmostEqual(d, 0)
        d = probability_to_distance(3.0 / 4.0)
        self.assertEqual(d, float('inf'))
    
    def test_conversion(self):
        probs = (0.0, 0.2, 0.4, 0.6)
        for p_expected in probs:
            d = probability_to_distance(p_expected)
            p_observed = distance_to_probability(d)
            self.assertAlmostEqual(p_expected, p_observed)

    def test_ml(self):
        sa = 'AAAAA'
        sbs = ('AAAAA', 'AAAAC', 'AAACC', 'AACCC')
        probs = (0.0, 0.2, 0.4, 0.6)
        for p_expected, sb in zip(probs, sbs):
            d_expected = probability_to_distance(p_expected)
            d_observed = get_ML_distance(sa, sb)
            self.assertAlmostEqual(d_expected, d_observed)

    def test_fisher_theano_f(self):
        mu = 0.6
        N = 4
        t = 1.4
        obj_plain = FisherInformation(mu, N)
        obj_theano = FisherInformationTheano(mu, N)
        self.assertAlmostEqual(obj_plain(t), obj_theano(t))

    def test_fisher_theano_df(self):
        mu = 0.6
        N = 4
        t = 1.4
        obj_plain = FisherInformation(mu, N)
        obj_theano = FisherInformationTheano(mu, N)
        self.assertAlmostEqual(obj_plain.deriv(t), obj_theano.deriv(t))


if __name__ == '__main__':
    unittest.main()
