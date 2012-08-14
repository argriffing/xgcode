"""
Compute fixation probability for a simple W-F model with selection.

Do this by solving a Markov chain numerically.
W-F is Wright-Fisher.
The fitness of B is 1+s, where s is determined by the user,
and the fitness of b is 1.
According to Haldane in 1927,
the fixation probability of a new advantageous allele is about 2s.
According to Kimura,
the fixation probability of an advantageous alelle is about
(1 - exp(-2*s*p*N)) / (1 - exp(-2*s*N))
where N is an effective population size,
1+s is the fitness of the advantageous allele,
and p is the initial frequency of the advantageous allele.
"""

from StringIO import StringIO
import math
from math import exp

import numpy as np
from scipy import linalg
from scipy import interpolate

import Form
import FormOut
import MatrixUtil
import StatsUtil

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('nB', 'number of B alleles', 1, low=0),
            Form.Integer('nb', 'number of b alleles', 400, low=0),
            Form.Float('s', 'positive when B is more fit', 0.01,
                low_exclusive=-1),
            ]

def get_form_out():
    return FormOut.Report()

def solve(npop, s):
    """
    @param npop: population
    @param s: selection
    @return: vector of probabilities
    """
    nstates = npop + 1
    fB = 1.0 + s
    fb = 1.0
    # compute the transition matrix 
    P = np.zeros((nstates, nstates))
    for i in range(nstates):
        nB_initial = i
        pB = (nB_initial * fB) / (nB_initial * fB + (npop - nB_initial) * fb)
        for j in range(nstates):
            nB_final = j
            log_p = StatsUtil.binomial_log_pmf(nB_final, npop, pB)
            P[i, j] = math.exp(log_p)
    # Put the puzzle into the form Ax=b
    # so that it can be solved by a generic linear solver.
    A = P - np.eye(nstates)
    b = np.zeros(nstates)
    # Adjust the matrix to disambiguate absorbing states.
    A[0, 0] = 1.0
    A[npop, npop] = 1.0
    b[0] = 0.0
    b[npop] = 1.0
    # Solve Ax=b for x.
    x = linalg.solve(A, b)
    # Return the probability of fixation of the B allele.
    return x

def get_response_content(fs):
    npop = fs.nB + fs.nb
    nstates = npop + 1
    # Check for minimum population size.
    if npop < 1:
        raise ValueError('there should be at least one individual')
    # Check the complexity;
    # solving a system of linear equations takes about n^3 effort.
    if nstates ** 3 > 1e8:
        raise ValueError('sorry this population size is too large')
    # Compute the exact probability of fixation of B.
    p_fixation = solve(npop, fs.s)[fs.nB]
    # Kimura approximation numerator and denominator.
    k_top = 1 - exp(-2*fs.s*fs.nB)
    k_bot = 1 - exp(-2*fs.s*npop)
    out = StringIO()
    print >> out, 'probability of eventual fixation (as opposed to extinction)'
    print >> out, 'of allele B in the population:'
    print >> out, p_fixation
    print >> out
    print >> out, 'Kimura would give the approximation'
    print >> out, k_top / k_bot
    print >> out
    if fs.nB == 1:
        print >> out, 'Haldane would give the approximation'
        print >> out, 2 * fs.s
        print >> out
    # Compute low-population approximations of probability of fixation of B.
    pB = fs.nB / float(fs.nB + fs.nb)
    for nsmall, name in (
            (10, 'low population size'),
            (20, 'medium population size'),
            ):
        if nsmall >= npop:
            continue
        s_small = fs.s * npop / float(nsmall)
        # Compute all low-population approximations.
        x = solve(nsmall, s_small)
        f_linear = interpolate.interp1d(range(nsmall+1), x, kind='linear')
        f_cubic = interpolate.interp1d(range(nsmall+1), x, kind='cubic')
        print >> out, 'linearly interpolated %s (N=%s)' % (name, nsmall)
        print >> out, 'approximation of probability of eventual'
        print >> out, 'fixation (as opposed to extinction)'
        print >> out, 'of allele B in the population:'
        print >> out, f_linear(pB*nsmall)
        print >> out
        print >> out, 'cubic interpolation:'
        print >> out, f_cubic(pB*nsmall)
        print >> out
    return out.getvalue()

