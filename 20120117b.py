"""
Investigate entrywise matrix operations related to Markov mutual information.

Computation of the decreasing mutual information as a function of time
can be written as either a bunch of nested summations
or by using matrix notation.
In the matrix notation, Hadamard products are required.
These products do not have as familiar properties as
standard matrix multiplication, especially with respect to eigendecomposition,
so I am doing some numerical computations to get a better feel
for what is going on.
The larger purpose is to get a better understanding of
asymptotics of mutual information through the spectrum
of a continuous time Markov rate matrix.
"""

from StringIO import StringIO
import math

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import divtime
import ctmcmi

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates', 'number of states', 4, low=2, high=9)]
            #Form.Float('divtime', 'arbitrary large-ish divergence time',
                #'3', low_exclusive=0)]
            #Form.Float('delta', 'vanishing time delta',
                #'0.0001', high_exclusive=0.25)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    out = StringIO()
    np.set_printoptions(linewidth=200)
    # get the user defined variables
    n = fs.nstates
    # sample a random rate matrix
    v = divtime.sample_distribution(n)
    S = divtime.sample_symmetric_rate_matrix(n)
    R = mrate.to_gtr_halpern_bruno(S, v)
    # get some properties of the rate matrix and its re-symmetrization
    S = mrate.symmetrized(R)
    distn = mrate.R_to_distn(R)
    w, U = np.linalg.eigh(S)
    D = np.diag(U.T[-1])**2
    D_inv = np.diag(np.reciprocal(U.T[-1]))**2
    for t in (1.0, 2.0):
        P = scipy.linalg.expm(R*t)
        M = ndot(D**.5, scipy.linalg.expm(S*t), D**.5)
        M_star = ndot(D_inv**.5, scipy.linalg.expm(S*t), D_inv**.5)
        M_star_log = np.log(M_star)
        M_star_log_w, M_star_log_U = np.linalg.eigh(M_star_log)
        E = M * np.log(M_star)
        E_w, E_U = np.linalg.eigh(E)
        print >> out, 't:'
        print >> out, t
        print >> out
        print >> out, 'randomly sampled rate matrix R'
        print >> out, R
        print >> out
        print >> out, 'symmetrized matrix S'
        print >> out, S
        print >> out
        print >> out, 'stationary distribution diagonal D'
        print >> out, D
        print >> out
        print >> out, 'R = D^-1/2 S D^1/2'
        print >> out, ndot(D_inv**.5, S, D**.5)
        print >> out
        print >> out, 'probability matrix e^(R*t) = P'
        print >> out, P
        print >> out
        print >> out, 'P = D^-1/2 e^(S*t) D^1/2'
        print >> out, ndot(D_inv**.5, scipy.linalg.expm(S*t), D**.5)
        print >> out
        print >> out, 'pairwise distribution matrix M'
        print >> out, 'M = D^1/2 e^(S*t) D^1/2'
        print >> out, M
        print >> out
        print >> out, 'sum of entries of M'
        print >> out, np.sum(M)
        print >> out
        print >> out, 'M_star = D^-1/2 e^(S*t) D^-1/2'
        print >> out, M_star
        print >> out
        print >> out, 'entrywise logarithm logij(M_star)'
        print >> out, np.log(M_star)
        print >> out
        print >> out, 'Hadamard product M o logij(M_star) = E'
        print >> out, E
        print >> out
        print >> out, 'spectrum of M:'
        print >> out, np.linalg.eigvalsh(M)
        print >> out
        print >> out, 'spectrum of logij(M_star):'
        print >> out, M_star_log_w
        print >> out
        print >> out, 'corresponding eigenvectors of logij(M_star) as columns:'
        print >> out, M_star_log_U
        print >> out
        print >> out, 'spectrum of E:'
        print >> out, E_w
        print >> out
        print >> out, 'corresponding eigenvectors of E as columns:'
        print >> out, E_U
        print >> out
        print >> out, 'entrywise square roots of stationary distribution:'
        print >> out, np.sqrt(v)
        print >> out
        print >> out, 'sum of entries of E:'
        print >> out, np.sum(E)
        print >> out
        print >> out, 'mutual information:'
        print >> out, ctmcmi.get_mutual_information(R, t)
        print >> out
        print >> out
    return out.getvalue()

