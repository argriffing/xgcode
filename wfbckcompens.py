"""
This is named in analogy with wffwdcompens.

The name of the module wffwdcompens means wright-fisher forward compensatory,
whereas the name of this module wfbckcompens means wright-fisher backward
compensatory in analogy.
The difference between the modules is that this module uses
linear algebra operations on matrices whereas the other module
uses forward simulation on state vectors.
"""

import numpy as np
from scipy import linalg

import MatrixUtil
import hittingtime

def get_substitution_info(P):
    """
    Get properties of an entire compensatory substitution.
    Unlike the type 1 and type 2 info functions,
    returns to AB fixation are counted towards the substitution path length.
    @param P: a huge transition matrix which is not modified
    @return: expectation and variance of compensatory substitution time
    """
    nstates = len(P)
    p = [i for i in range(nstates) if i != 3]
    Q = P[p, :][:, p]
    n = nstates - 1
    t = linalg.solve(np.eye(n) - Q, np.ones(n))
    v = 2*linalg.solve(np.eye(n) - Q, t) - t * (t + 1)
    return t[0], v[0]

def get_type_2_info(P):
    """
    The expected time for a type 2 event is computed as follows.
    It is the expected number of steps from AB to ab
    conditional on not entering the states AB, Ab, or aB.
    It should also include a bit of exponential delay that it takes
    to leave the final fixed AB state before embark.
    @param P: a huge transition matrix which is not modified
    @return: expectation and variance of compensatory substitution time
    """
    MatrixUtil.assert_transition_matrix(P)
    nstates = len(P)
    # define index sequences
    plain = range(4, nstates)
    forbidden = [0, 1, 2]
    target = [3]
    #
    H = hittingtime.get_conditional_transition_matrix(
            P, plain, forbidden, target)
    t = hittingtime.get_absorption_time(
            H, plain+forbidden, target)
    v = hittingtime.get_absorption_variance(
            H, plain+forbidden, target)
    #
    t0 = t[0]
    v0 = v[0]
    # add a geometric rv that depends on probability of leaving fixed AB
    p = 1 - P[0, 0]
    t0 += (1 - p) / p
    v0 += (1 - p) / (p*p)
    #
    return t0, v0

def get_type_1_info(P):
    """
    The expected time for a type 1 event is computed as follows.
    It is the sum of two expected times.
    The first time is the expected number of steps from AB to either Ab or aB
    conditional on not entering the states AB or ab.
    The second time is the expected number of steps from Ab to ab
    conditional on not entering AB.
    Note that this formulation depends on the assumption
    that the process associated with the first
    step is equally likely to end up in Ab as in aB.
    It should also include a bit of exponential delay that it takes
    to leave the final fixed AB state before embark.
    @param P: a huge transition matrix which is not modified
    @return: expectation and variance of compensatory substitution time
    """
    MatrixUtil.assert_transition_matrix(P)
    nstates = len(P)
    # get the expected time for the first stage
    plain = range(4, nstates)
    forbidden = [0, 3]
    target = [1, 2]
    H = hittingtime.get_conditional_transition_matrix(
            P, plain, forbidden, target)
    t = hittingtime.get_absorption_time(
            H, plain+forbidden, target)
    v = hittingtime.get_absorption_variance(
            H, plain+forbidden, target)
    t_first = t[0]
    v_first = v[0]
    # get the expected time for the second stage
    plain = [1, 2] + range(4, nstates)
    forbidden = [0]
    target = [3]
    H = hittingtime.get_conditional_transition_matrix(
            P, plain, forbidden, target)
    t = hittingtime.get_absorption_time(
            H, plain+forbidden, target)
    v = hittingtime.get_absorption_variance(
            H, plain+forbidden, target)
    t_second = t[1]
    v_second = v[1]
    # add a geometric rv that depends on probability of leaving fixed AB
    p = 1 - P[0, 0]
    t_third = (1 - p) / p
    v_third = (1 - p) / (p*p)
    # return the moments of the distribution accounting for both stages
    return t_first + t_second + t_third, v_first + v_second + v_third

