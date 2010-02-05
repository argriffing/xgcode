"""Functions dealing with the F84 nucleotide rate matrix.

The F84 evolutionary model is defined in the paper
"Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods"
by Ziheng Yang in J Mol Evol 1994.
"""

from StringIO import StringIO
import unittest
import math

import scipy.optimize

import MatrixUtil
import RateMatrix
import Fasta
import PairLikelihood


g_transitions = (('C', 'T'), ('T', 'C'), ('A', 'G'), ('G', 'A'))

def create_rate_matrix(kappa, nt_distribution):
    """
    @param kappa: adjusts for the transition rate differing from the transversion rate
    @param nt_distribution: ordered ACGT nucleotide probabilities
    @return: a rate matrix object with one expected nucleotide substitution per time unit
    """
    # make some assertions about the distribution
    for p in nt_distribution:
        assert p >= 0
    assert len(nt_distribution) == 4
    assert RateMatrix.almost_equal(sum(nt_distribution), 1.0)
    # define some intermediate variables
    A, C, G, T = nt_distribution
    R = float(A + G)
    Y = float(C + T)
    # make some more assertions about the distribution and about kappa
    assert A+G > 0
    assert C+T > 0
    assert kappa > max(-Y, -R)
    # get the normalization constant
    normalization_constant = 4*T*C*(1 + kappa/Y) + 4*A*G*(1 + kappa/R) + 4*Y*R
    # adjust the normalization constant to correct what might be an error in the paper
    normalization_constant /= 2
    # define the dictionary rate matrix
    dict_rate_matrix = {}
    for source_index, source in enumerate('ACGT'):
        for sink_index, sink in enumerate('ACGT'):
            key = (source, sink)
            coefficient = 1.0
            if key in g_transitions:
                coefficient = 1 + kappa / (nt_distribution[source_index] + nt_distribution[sink_index])
            dict_rate_matrix[key] = coefficient * nt_distribution[sink_index] / normalization_constant
    for source in 'ACGT':
        dict_rate_matrix[(source, source)] = -sum(dict_rate_matrix[(source, sink)] for sink in 'ACGT' if source != sink)
    # convert the dictionary rate matrix to a row major rate matrix
    row_major = MatrixUtil.dict_to_row_major(dict_rate_matrix, 'ACGT', 'ACGT')
    # return the rate matrix object
    rate_matrix_object = RateMatrix.RateMatrix(row_major, 'ACGT')
    expected_rate = rate_matrix_object.get_expected_rate()
    if not RateMatrix.almost_equal(expected_rate, 1.0):
        assert False, 'the rate is %f but should be 1.0' % expected_rate
    return rate_matrix_object

def get_closed_form_estimates(sequence_pair):
    """
    This will fail if some of the nucleotides are unobserved.
    This formula is from Ziheng Yang JME 1994.
    @param sequence_pair: a pair of upper case nucleotide sequences
    @return: (distance, kappa, A, C, G, T)
    """
    # assert that the input is reasonable
    assert len(sequence_pair) == 2
    # convert the generic sequences to strings
    sa = ''.join(list(sequence_pair[0])).upper()
    sb = ''.join(list(sequence_pair[1])).upper()
    st = sa+sb
    # assert that the sequences are long enough and are the same length
    assert len(sa) > 1
    assert len(sa) == len(sb)
    # assert that the strings are composed of valid nucleotides
    assert set(st) <= set('ACGT')
    # get the estimated nucleotide distribution
    n = len(sa)
    assert len(st) == 2*n
    A = st.count('A') / float(2*n)
    C = st.count('C') / float(2*n)
    G = st.count('G') / float(2*n)
    T = st.count('T') / float(2*n)
    # get some intermediate variables
    R = A + G
    Y = C + T
    if not R or not Y:
        assert False, 'expected some purines and some pyrimidines'
    # get the fraction of sites with transitional and with transversional differences
    P = 0
    Q = 0
    for a, b in zip(sa, sb):
        key = (a, b)
        if a == b:
            continue
        elif key in g_transitions:
            P += 1
        else:
            Q += 1
    P /= float(n)
    Q /= float(n)
    # get some more intermediate variables
    A_num = 2*(T*C + A*G) + 2*(T*C*R/Y + A*G*Y/R) * (1 - Q/(2*Y*R)) - P
    A_den = 2*(T*C/Y + A*G/R)
    A_total = A_num / A_den
    B_total = 1 - Q/(2*Y*R)
    if A_total <= 0 or B_total <= 0:
        assert False, 'range error in an intermediate calculation'
    alpha = -math.log(A_total)/2
    beta = -math.log(B_total)/2
    # estimate kappa
    kappa = alpha / beta - 1
    # estimate the distance
    distance = beta * (4*T*C*(1+kappa/Y) + 4*A*G*(1+kappa/R) + 4*Y*R)
    # return the estimates
    return (distance, kappa, A, C, G, T)


def parameters_to_distribution(parameters):
    """
    @param parameters: these N-1 degrees of freedom define a distribution over N states
    """
    den = 1.0 + sum(math.exp(p) for p in parameters)
    dist = [1.0 / den]
    for p in parameters:
        dist.append(math.exp(p) / den)
    return dist

class Objective:
    """
    This class defines a function object used for testing maximum likelihood estimation.
    This is currently used only for testing.
    An instance will be called using scipy.optimize.fmin(...).
    """
    def __init__(self, sequence_pair):
        """
        @param sequence_pair: the observed data
        """
        self.sequence_pair = sequence_pair

    def get_initial_parameters(self):
        """
        @return: a vector of parameters
        """
        # get an initial guess
        distance, kappa, A, C, G, T = get_closed_form_estimates(self.sequence_pair)
        wC, wG, wT = (0, 0, 0)
        # perturb the initial guess just to show that it is arbitrary
        perturbation = .01
        theta = distance + perturbation, kappa + perturbation, wC, wG, wT
        # return the initial guess
        return theta

    def __call__(self, theta):
        """
        @param theta: the vector of estimated parameters
        @return: the negative log likelihood to be minimized
        """
        # unpack the parameters
        distance, kappa, wC, wG, wT = theta
        nt_distribution = parameters_to_distribution((wC, wG, wT))
        # make the rate matrix
        model = create_rate_matrix(kappa, nt_distribution)
        # get the likelihood
        log_likelihood = PairLikelihood.get_log_likelihood(distance, self.sequence_pair, model)
        return -log_likelihood


class TestF84(unittest.TestCase):
    
    def testMLE(self):
        """
        Assert that there is no syntax error in the MLE estimation code.
        """
        sa = 'AAACG'
        sb = 'TTACG'
        pair = (sa, sb)
        distance, kappa, A, C, G, T = get_closed_form_estimates(pair)

    def testCreate(self):
        """
        Assert that there is no syntax error in the matrix creation code.
        """
        nt_distribution = (.25, .25, .25, .25)
        kappa = 1.0
        rate_matrix_object = create_rate_matrix(kappa, nt_distribution)

    def testNormalization(self):
        """
        Assert that the matrix is properly normalized.
        """
        nt_distribution = (.1, .2, .3, .4)
        kappa = 2.0
        rate_matrix_object = create_rate_matrix(kappa, nt_distribution)
        self.assertAlmostEqual(1.0, rate_matrix_object.get_expected_rate())

    def testLikelihood(self):
        """
        Assert that no errors occur during the analysis
        """
        # define a simple (but not completely degenerate) alignment
        sa = 'AAAACCCCGGGGTTAA'
        sb = 'GAAACCTCGGCGTAAA'
        sequence_pair = (sa, sb)
        # get estimates according to an analytical formula which is not necessarily the mle
        distance_mle, kappa_mle, A_mle, C_mle, G_mle, T_mle = get_closed_form_estimates((sa, sb))
        nt_distribution_mle = (A_mle, C_mle, G_mle, T_mle)
        rate_matrix_object = create_rate_matrix(kappa_mle, nt_distribution_mle)
        log_likelihood_mle = PairLikelihood.get_log_likelihood(distance_mle, sequence_pair, rate_matrix_object)
        # get the maximum likelihood estimates according to a numeric optimizer.
        f = Objective((sa, sb))
        values = list(f.get_initial_parameters())
        result = scipy.optimize.fmin(f, values, ftol=.0000000001, disp=0)
        distance_opt, kappa_opt, wC_opt, wG_opt, wT_opt = result
        nt_distribution_opt = parameters_to_distribution((wC_opt, wG_opt, wT_opt))
        rate_matrix_object = create_rate_matrix(kappa_opt, nt_distribution_opt)
        log_likelihood_opt = PairLikelihood.get_log_likelihood(distance_opt, sequence_pair, rate_matrix_object)


if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestF84)
        unittest.TextTestRunner(verbosity=2).run(suite)

