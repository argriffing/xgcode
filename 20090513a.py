"""Check an example of partial autocorrelation.

Say that there is a discrete time process where the deltas are alternating in sign.
What does this imply about the lag 1 and lag 2 partial autocorrelations?
I am doing this analysis out of curiosity based on a poster by Reed Cartwright.
"""

import StringIO
import random

from SnippetUtil import HandlingError
import Form

def get_form():
    """
    @return: the body of a form
    """
    return [Form.Integer('n', 'sequence length', 100, low=3, high=1000)]

def get_pairs(seq, lag):
    """
    @param seq: a sequence of numbers
    @param lag: the difference between paired indices
    @return: pairs of elements of the sequence, separated by the given lag
    """
    return zip(seq[:-lag], seq[lag:])

def get_autocorrelation(seq, lag):
    """
    @param seq: a sequence of numbers
    @param lag: the difference between paired indices
    @return: a floating point autocorrelation between -1 and 1
    """
    n = len(seq)
    mu = sum(seq) / float(n)
    numerator = sum((a-mu)*(b-mu) for a, b in get_pairs(seq, lag))
    denominator = sum((a-mu)*(a-mu) for a in seq)
    return numerator / denominator

def get_pacf_1(seq):
    """
    @param seq: a sequence of numbers
    @return: the first partial autocorrelation of the sequence
    """
    return get_autocorrelation(seq, 1)

def get_pacf_2(seq):
    """
    @param seq: a sequence of numbers
    @return: the second partial autocorrelation of the sequence
    """
    rho_1 = get_autocorrelation(seq, 1)
    rho_2 = get_autocorrelation(seq, 2)
    pacf_2 = (rho_2 - rho_1**2) / (1 - rho_1**2)
    return pacf_2

def sample_sequence_same_sign(n):
    """
    @param n: the length of the sequence
    @return: a sequence of n numbers
    """
    seq = []
    current = 0
    for i in range(100):
        delta = random.randrange(1, 3)
        sign = 1
        current += delta * sign
        seq.append(current)
    return seq

def sample_sequence_alternating_sign(n):
    """
    @param n: the length of the sequence
    @return: a sequence of n numbers
    """
    seq = []
    current = 0
    for i in range(100):
        delta = random.randrange(1, 3)
        sign = 1 if (i % 2) else -1
        current += delta * sign
        seq.append(current)
    return seq

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # define a test sequence from page 22 of
    # "Time series analysis: univariate and multivariate methods"
    # by Wei
    test_seq = [13, 8, 15, 4, 4, 12, 11, 7, 14, 12]
    assert len(test_seq) == 10
    # define the response
    out = StringIO.StringIO()
    # show test results
    print >> out, 'test results:'
    print >> out, 'expected autocorrelation of lag 1: -0.188'
    print >> out, 'observed autocorrelation of lag 1:', get_autocorrelation(test_seq, 1)
    print >> out, 'expected autocorrelation of lag 2: -0.201'
    print >> out, 'observed autocorrelation of lag 2:', get_autocorrelation(test_seq, 2)
    print >> out, 'expected partial autocorrelation of lag 1: -0.188'
    print >> out, 'observed partial autocorrelation of lag 1:', get_pacf_1(test_seq)
    print >> out, 'expected partial autocorrelation of lag 2: -0.245'
    print >> out, 'observed partial autocorrelation of lag 2:', get_pacf_2(test_seq)
    print >> out
    # show simulation results
    seq = sample_sequence_alternating_sign(fs.n)
    print >> out, 'simulation results for deltas with alternating sign:'
    print >> out, 'observed autocorrelation of lag 1:', get_autocorrelation(seq, 1)
    print >> out, 'observed autocorrelation of lag 2:', get_autocorrelation(seq, 2)
    print >> out, 'observed partial autocorrelation of lag 1:', get_pacf_1(seq)
    print >> out, 'observed partial autocorrelation of lag 2:', get_pacf_2(seq)
    print >> out
    # show simulation results
    seq = sample_sequence_same_sign(fs.n)
    print >> out, 'simulation results for deltas with same sign:'
    print >> out, 'observed autocorrelation of lag 1:', get_autocorrelation(seq, 1)
    print >> out, 'observed autocorrelation of lag 2:', get_autocorrelation(seq, 2)
    print >> out, 'observed partial autocorrelation of lag 1:', get_pacf_1(seq)
    print >> out, 'observed partial autocorrelation of lag 2:', get_pacf_2(seq)
    print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
