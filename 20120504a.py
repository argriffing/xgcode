r"""
Compute the number of significant principal components.

This uses a variant of Tracy-Widom theory from Patterson et al.
Let \( M_{m \times n} \) be a matrix where \( m &lt n \).
Let \( X = M M^T \) be the matrix whose \( m \) eigenvalues we care about.
Then we want to find the number of 'significant' principal components
given \(m, n\) and the eigenvalues of \( X \).
For the purposes of this script,
if the number of eigenvalues given is less than \( m \)
then the rest of the eigenvalues are assumed to be near zero.
The default values are taken from Table 3 of Patterson et al.
except presumably they had more than the ten nonzero eigenvalues
printed in their table.
For the Tracy-Widom cumulative density function table,
see 30983-cdf-for-tracy-widom-tw1-distribution
at the matlab file exchange.
"""

from StringIO import StringIO
import math
from math import sqrt

import numpy as np
import scipy
from scipy import interpolate
from scipy import optimize

import Form
import FormOut
import Util
import const

g_tw1_csv = const.read('20120505a')

g_m = 189

g_n = 2790

g_eigenvalues = """
22.36 
8.20 
5.09 
3.81 
3.33 
2.09 
1.89 
1.44 
1.30 
1.27 
""".strip().splitlines()

# According to Patterson et al.
# the critical points of the Tracy-Widom density
# corresponding to some conventional significance levels are
# (p=0.05, x=0.9794),
# (p=0.01, x=2.0236),
# (p=0.001, x=3.2730).

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('m', 'm', g_m, low=2),
            Form.Integer('n', 'n', g_n, low=2),
            Form.Sequence('eigenvalues', 'eigenvalues', g_eigenvalues),
            Form.Float('significance_level', 'significance level',
                '0.05', low_exclusive=0, high_exclusive=1),
            Form.RadioGroup('method', 'method from Patterson et al.', [
                Form.RadioItem('method_eq7', 'Eq. 7'),
                Form.RadioItem('method_eq8', 'Eq. 8'),
                Form.RadioItem('method_a_test_for_pop',
                    'the section "A Test for Population Structure"', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

class Err:
    """
    The root of this callable is the critical value.
    """
    def __init__(self, cdf, sig_level):
        self.cdf = cdf
        self.sig_level = sig_level
    def __call__(self, x):
        return self.cdf(x) - (1 - self.sig_level)

def get_critical_value(significance_level):
    lines = Util.get_stripped_lines(g_tw1_csv.splitlines())
    s_pairs = [line.split(',') for line in lines]
    xy_pairs = [(float(x), float(y)) for x, y in s_pairs]
    X, Y = zip(*xy_pairs)
    cdf = scipy.interpolate.UnivariateSpline(X, Y, s=0)
    error_function = Err(cdf, significance_level)
    return scipy.optimize.brentq(error_function, X[0], X[-1])

def _tw(m, n, lambda_like):
    """
    @param m: rows of the data matrix
    @param n: columns of the data matrix
    @param lambda_like: like a principal eigenvalue
    @return: an approximation of the Tracy-Widom statistic
    """
    alpha = sqrt(n-1) + sqrt(m)
    mu = (alpha * alpha) / n
    sigma = (alpha / n) * (1/sqrt(n-1) + 1/sqrt(m))**(1./3.)
    return (lambda_like - mu) / sigma

def get_eq_7_statistic(m, n, W):
    """
    @param m: number of rows of the data matrix
    @param n: number of columns of the data matrix
    @param lambda_1: first principal eigenvalue of the wishart-like matrix
    @return: an approximation of the Tracy-Widom statistic
    """
    lambda_1 = W[0]
    return _tw(m, n, lambda_1)

def get_eq_8_statistic(m, n, W):
    """
    @param m: number of rows of the data matrix
    @param n: number of columns of the data matrix
    @param W: decreasing eigenvalues of the wishart-like matrix
    @return: an approximation of the Tracy-Widom statistic
    """
    L1 = m * (W[0] / np.sum(W))
    return _tw(m, n, L1)

def get_effective_n(m, W):
    """
    Equation 10 from Patterson et al.
    The un-numbered equation that follows equation 10 gives an estimated
    variance of the entries of M, but this is not used.
    @param m: number of rows of the data matrix
    @param W: decreasing eigenvalues of the wishart-like matrix
    @return: an effective number of columns, and an estimated sigma
    """
    n_prime_numerator = (m + 1) * np.sum(W)**2
    n_prime_denominator = ((m - 1) * np.dot(W, W)) - np.sum(W)**2
    n_prime = n_prime_numerator / n_prime_denominator
    return n_prime

def get_a_test_for_population_structure_statistic(m, n, W):
    """
    @param m: number of rows of the data matrix
    @param n: number of columns of the data matrix
    @param W: decreasing eigenvalues of the wishart-like matrix
    @return: an approximation of the Tracy-Widom statistic
    """
    m_prime = m - 1
    n_prime = get_effective_n(m, W)
    # if the method of moments estimate of the effective number of columns
    # is greater than the actual number of columns,
    # then use the actual number of columns.
    n_prime = min(n_prime, n)
    ell = (m_prime * W[0]) / np.sum(W)
    # FIXME should the following line use m_prime instead of m?
    return _tw(m, n_prime, ell)

def get_response_content(fs):
    # get the user data
    m = fs.m
    n = fs.n
    if not fs.eigenvalues:
        raise ValueError('expected at least one eigenvalue')
    W = [float(x) for x in fs.eigenvalues]
    W = np.array(sorted(W, reverse=True))
    significance_level = fs.significance_level
    critical_value = get_critical_value(significance_level)
    if m < len(W):
        raise ValueError('the number of eigenvalues should be at most m')
    # define the function to compute the Tracy-Widom statistic
    if fs.method_eq7:
        f = get_eq_7_statistic
    elif fs.method_eq8:
        f = get_eq_8_statistic
    elif fs.method_a_test_for_pop:
        f = get_a_test_for_population_structure_statistic
    # compute some stuff using the Patterson et al.
    x_values = []
    nsignificant = 0
    for k, eigenvalue in enumerate(W):
        m_remaining = m - k
        x = f(m_remaining, n, W[k:])
        x_values.append(x)
        if x < critical_value:
            break
        else:
            nsignificant += 1
    # begin the output
    out = StringIO()
    # show the number of significant principal components
    print >> out, nsignificant, 'significant principal components'
    print >> out
    print >> out, 'specified significance level:', significance_level
    print >> out, 'corresponding critical value:', critical_value
    print >> out
    # show the eigenvalues and statistics
    print >> out, 'number\teigenvalue\tstatistic'
    for k, (eigenvalue, x) in enumerate(zip(W, x_values)):
        print >> out, '%s\t%s\t%s' % (k+1, eigenvalue, x)
    return out.getvalue()

