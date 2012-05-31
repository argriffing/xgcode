r"""
Look for selection that maximizes Fisher information for a given divergence.

The mutation process is assumed to be a
site-independent 3-site 2-state-per-site process
with uniform stationary distribution.
The selection is assumed to work according to a formula used
by Sella and Hirsh and Halpern and Bruno.
"""

from StringIO import StringIO
import argparse
import math
import time
import random
from itertools import product

import numpy as np
import scipy
from scipy import linalg
from scipy import optimize

import Form
import FormOut
import ctmcmi
import mrate
import divtime
import cheeger
import MatrixUtil
from MatrixUtil import ndot

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('t', 'divergence time', '2.0', low_exclusive=0),
            Form.RadioGroup('infotype', 'information type', [
                Form.RadioItem('infotype_fi', 'Fisher information', True),
                Form.RadioItem('infotype_mi', 'mutual information')]),
            ]

def get_form_out():
    return FormOut.Report('results')

def get_site_independent_distn(v, nsites):
    n = len(v)
    v_out = np.zeros(n**nsites)
    for values in product(range(n), repeat=nsites):
        index = sum(k * (n**i) for i, k in enumerate(values))
        p = 1
        for value in values:
            p *= v[value]
        v_out[index] = p
    return v_out

def get_site_independent_process(R, nsites):
    """
    @param R: single site process
    @param nsites: return a rate matrix that acts on this many sites
    @return: a site independent rate matrix
    """
    n = len(R)
    nprod = n**nsites
    Q = np.zeros((nprod, nprod))
    for a_values in product(range(n), repeat=nsites):
        a_index = sum(k * (n**i) for i, k in enumerate(a_values))
        for b_values in product(range(n), repeat=nsites):
            b_index = sum(k * (n**i) for i, k in enumerate(b_values))
            diffs = [(i, a, b) for i, (a, b) in enumerate(zip(
                a_values, b_values)) if a != b]
            if len(diffs) == 1:
                i, a, b = diffs[0]
                Q[a_index, b_index] = R[a, b]
    Q -= np.diag(np.sum(Q, axis=1))
    return Q

#FIXME this is copypasted from another web script
class OptDep:
    """
    This is for numerical optimization.
    """
    def __init__(self, R, v, t, f_info, f_selection):
        """
        @param R: two site independent mutation rate matrix
        @param v: stationary distribution of R
        @param t: time
        @param f_info: info function that takes a rate matrix and a time
        @param f_selection: takes a rate matrix and a current and target distn
        """
        self.R = R
        self.v = v
        self.t = t
        self.f_info = f_info
        self.f_selection = f_selection
    def get_process(self, X):
        """
        @param X: some log ratio probabilities
        @return: the mut-sel balance rate matrix and stationary distn
        """
        if len(X) != len(self.R) - 1:
            raise ValueError('state space size mismatch')
        # define the mutation and balance energies
        u = -(np.log(self.v) - math.log(self.v[0]))
        v = np.array([0.0] + X.tolist())
        # apply the selection using the stable halpern bruno function
        Q = self.f_selection(self.R, u, v)
        # get the stationary distribution for downstream
        qdistn = np.exp(-v)
        qdistn /= np.sum(qdistn)
        if np.any(np.isnan(Q)):
            print self.R
            print qdistn
            print Q
            raise ValueError('the rate matrix has nans')
        return Q, qdistn
    def __call__(self, X):
        """
        @param X: some log ratio probabilities
        @return: neg info value for minimization
        """
        Q, v_target = self.get_process(X)
        return -self.f_info(Q, v_target, self.t)

#FIXME this is copypasted from another web script
class OptDepNearParity(OptDep):
    """
    This assumes a very hardcoded mutation process.
    """
    def get_process(self, X):
        """
        @param X: unpacks to a single energy
        @return: the mut-sel balance rate matrix and stationary distn
        """
        g, = X
        # define the mutation and balance energies
        u = -(np.log(self.v) - math.log(self.v[0]))
        v = np.array([g, 0, 0, g, 0, g, g, 0])
        # apply the selection using the stable halpern bruno function
        Q = self.f_selection(self.R, u, v)
        # get the stationary distribution for downstream
        qdistn = np.exp(-v)
        qdistn /= np.sum(qdistn)
        if np.any(np.isnan(Q)):
            print self.R
            print qdistn
            print Q
            raise ValueError('the rate matrix has nans')
        return Q, qdistn

def do_search_opt_dep(opt_dep, df, nseconds):
    """
    @param opt_dep: a function object for numerical optimization
    @param df: the number of degrees of freedom
    @param nseconds: spend roughly this many seconds searching
    @return: a freeform multiline string report
    """
    # utility vars
    n = len(M)
    uniform_distn = np.ones(n, dtype=float) / n
    f_selection = mrate.to_gtr_hb_known_energies
    # init the search
    best_info = None
    best_state = None
    t0 = time.time()
    niter = 0
    while time.time() - t0 < nseconds:
        niter += 1
        # get the site-dependent mutation selection balance information
        X0 = np.random.randn(df)
        xopt = scipy.optimize.fmin(
                opt_dep, X0, maxiter=10000, maxfun=10000)
        info = -opt_dep(xopt)
        if best_info is None or info > best_info:
            Q_bal, v_bal = opt_dep.get_process(xopt)
            state = [
                    info,
                    v_bal,
                    Q_bal,
                    ]
            best_info = info
            best_state = state
    out = StringIO()
    print >> out, 'numerical maximization runs:', niter
    for thing in state:
        print >> out, str(thing)
    return out.getvalue().rstrip()

#TODO this is way too copypasted and hardcoded
def do_search_near_parity(M, t, f_info, nseconds):
    """
    @param M: mutation process assumed uniform stationary distribution
    @param t: the divergence time for the process
    @param f_info: information function
    @param nseconds: number of seconds allowed for search
    @return: search result text
    """
    # utility vars
    n = len(M)
    uniform_distn = np.ones(n, dtype=float) / n
    f_selection = mrate.to_gtr_hb_known_energies
    # init the search
    best_info = None
    best_state = None
    t0 = time.time()
    niter = 0
    while time.time() - t0 < nseconds:
        niter += 1
        # get the site-dependent mutation selection balance information
        dep_balance = OptDepNearParity(M, uniform_distn, t, f_info, f_selection)
        X0 = np.random.randn(1)
        xopt = scipy.optimize.fmin(
                dep_balance, X0, maxiter=10000, maxfun=10000)
        info = -dep_balance(xopt)
        if best_info is None or info > best_info:
            Q_bal, v_bal = dep_balance.get_process(xopt)
            state = [
                    info,
                    v_bal,
                    Q_bal,
                    ]
            best_info = info
            best_state = state
    out = StringIO()
    print >> out, 'numerical maximization runs:', niter
    for thing in state:
        print >> out, str(thing)
    return out.getvalue().rstrip()

def do_search(M, t, f_info, nseconds):
    """
    @param M: mutation process assumed uniform stationary distribution
    @param t: the divergence time for the process
    @param f_info: information function
    @param nseconds: number of seconds allowed for search
    @return: search result text
    """
    # utility vars
    n = len(M)
    uniform_distn = np.ones(n, dtype=float) / n
    f_selection = mrate.to_gtr_hb_known_energies
    # init the search
    best_info = None
    best_state = None
    t0 = time.time()
    niter = 0
    while time.time() - t0 < nseconds:
        niter += 1
        # get the site-dependent mutation selection balance information
        dep_balance = OptDep(M, uniform_distn, t, f_info, f_selection)
        X0 = np.random.randn(n - 1)
        xopt = scipy.optimize.fmin(
                dep_balance, X0, maxiter=10000, maxfun=10000)
        info = -dep_balance(xopt)
        if best_info is None or info > best_info:
            Q_bal, v_bal = dep_balance.get_process(xopt)
            state = [
                    info,
                    v_bal,
                    Q_bal,
                    ]
            best_info = info
            best_state = state
    out = StringIO()
    print >> out, 'iterations:', niter
    for thing in state:
        print >> out, str(thing)
    return out.getvalue().rstrip()


#TODO recode this with combobreaker
def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    # hardcode some stuff
    nseconds_per_constraint = 2
    nsites = 3
    f_selection = mrate.to_gtr_hb_known_energies
    # validate and store user input
    t = fs.t
    if fs.infotype_fi:
        f_info = divtime.get_fisher_info_known_distn_fast
    elif fs.infotype_mi:
        f_info = ctmcmi.get_mutual_info_known_distn
    # define the single site and multi site mutation rate matrix
    M_two_state = 0.5 * np.array([[-1, 1], [1, -1]])
    M_box = get_site_independent_process(M_two_state, nsites)
    # Define some sub-processes for which some
    # states have been removed by lethal selection.
    M_coil_in_the_box = 0.5 * np.array([
        [-2, 1, 0, 0, 0, 1],
        [1, -2, 1, 0, 0, 0],
        [0, 1, -2, 1, 0, 0],
        [0, 0, 1, -2, 1, 0],
        [0, 0, 0, 1, -2, 1],
        [1, 0, 0, 0, 1, -2]])
    M_snake_in_the_box = 0.5 * np.array([
        [-1, 1, 0, 0, 0],
        [1, -2, 1, 0, 0],
        [0, 1, -2, 1, 0],
        [0, 0, 1, -2, 1],
        [0, 0, 0, 1, -1]])
    # do the searches
    #uniform_distn = np.ones(n, dtype=float) / n
    out = StringIO()
    for M, name in (
            (M_box, 'cube near parity'),
            ):
        print >> out, 'selection constraint:', name
        print >> out, M
        print >> out, do_search_near_parity(
                M, t, f_info, nseconds_per_constraint)
        print >> out
    for M, name in (
            (M_box, 'cube'),
            (M_snake_in_the_box, 'snake-in-the-box'),
            (M_coil_in_the_box, 'coil-in-the-box')):
        print >> out, 'selection constraint:', name
        print >> out, M
        print >> out, do_search(
                M, t, f_info, nseconds_per_constraint)
        print >> out
    return out.getvalue()

