"""
Compute fixation probability for a simple W-F model with selection.

Do this by solving a Markov chain numerically.
W-F is Wright-Fisher.
It assumes random mating and selection acting on the children.
The relative fitness of BB is 1
and the relative fitness of bb is 1-s where s is determined by the user.
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
import wfengine

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

def solve(N_haploid, s):
    """
    @param N_haploid: population
    @param s: selection
    @return: vector of probabilities
    """
    if N_haploid % 2:
        raise ValueError('expected an even haploid population')
    N_diploid = N_haploid / 2
    #fB = 1.0 + s
    #fb = 1.0
    # compute the transition matrix 
    P = np.exp(wfengine.create_genic_diallelic(N_diploid, s*2))
    # define some boundary conditions
    #P[0, 0] = 0
    #P[0, 1] = 1
    #P[N_diploid, N_diploid] = 0
    #P[N_diploid, 1] = 1
    #
    # Put the puzzle into the form Ax=b
    # so that it can be solved by a generic linear solver.
    A = P - np.eye(N_haploid + 1)
    b = np.zeros(N_haploid + 1)
    # Adjust the matrix to disambiguate absorbing states.
    A[0, 0] = 1
    A[N_haploid, N_haploid] = 1
    b[0] = 0
    b[N_haploid] = 1
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

