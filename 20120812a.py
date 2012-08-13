"""
Compute fixation probability for a very simple Wright-Fisher model.

Do this by solving a Markov chain numerically.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import linalg

import Form
import FormOut
import MatrixUtil
import StatsUtil

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('nB', 'number of B alleles', 2, low=0),
            Form.Integer('nb', 'number of b alleles', 6, low=0),
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    npop = fs.nB + fs.nb
    nstates = npop + 1
    # Check the complexity;
    # solving a system of linear equations takes about n^3 effort.
    if nstates ** 3 > 1e6:
        raise ValueError('sorry this population size is too large')
    # Compute the transition matrix.
    # This assumes no mutation or selection or recombination.
    # It is pure binomial.
    P = np.zeros((nstates, nstates))
    for i in range(nstates):
        nB_initial = i
        for j in range(nstates):
            nB_final = j
            log_p = StatsUtil.binomial_log_pmf(
                    nB_final, npop, nB_initial / float(npop))
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
    # Print the solution.
    out = StringIO()
    print >> out, 'probability of eventual fixation (as opposed to extinction)'
    print >> out, 'of allele B in the population:'
    print >> out, x[fs.nB]
    return out.getvalue()

