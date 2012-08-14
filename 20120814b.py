"""
Compute fixation probability for a diploid W-F model with selection.

Do this by solving a Markov chain numerically.
W-F is Wright-Fisher.
The homozygous aa type has fitness 1 by convention.
The other two types have fitness 1+s_type where the types are AA and Aa.
A continuous approximation is provided by Kimura in
"On the Probability of Fixation of Mutant Genes in a Population."
"""

from StringIO import StringIO
import math
from math import exp

import numpy as np
from scipy import linalg
from scipy import interpolate
from scipy import integrate

import Form
import FormOut
import MatrixUtil
import StatsUtil
import Util

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('nAA', 'number of AA genotypes', 0, low=0),
            Form.Integer('nAa', 'number of Aa genotypes', 1, low=0),
            Form.Integer('naa', 'number of aa genotypes', 20, low=0),
            Form.Float('sAA', 'AA selection value', 0.02, low_exclusive=-1),
            Form.Float('sAa', 'Aa selection value', 0.01, low_exclusive=-1),
            ]

def get_form_out():
    return FormOut.Report()

def get_state_space_size(npop):
    """
    How many ways to distribute N individuals among k groups.
    This is another way to check that the state space is not incorrect.
    See example 4.15 of i.stanford.edu/~ullman/focs/ch04.pdf
    """
    N = npop
    k = 3
    return Util.choose(N + k - 1, k - 1)

def gen_population_compositions(npop):
    nstates = npop + 1
    for i in range(nstates):
        for j in range(nstates - i):
            yield (i, j, npop - (i+j))

def get_child_distn(i, j):
    """
    The indices i and j are the two parent genotypes.
    Return a distribution over the child genotypes.
    Aggregate variants are 0:AA, 1:Aa, 2:aa
    @param i: aggregate variant index
    @param j: aggregate variant index
    @return: two aggregate variant indices
    """
    index_to_string = ('AA', 'Aa', 'aa')
    string_to_index = {
            'AA' : 0,
            'Aa' : 1,
            'aA' : 1,
            'aa' : 2,
            }
    si, sj = index_to_string[i], index_to_string[j]
    child_distn = np.zeros(3)
    for i in range(2):
        for j in range(2):
            s_child = si[i] + sj[j]
            child_distn[string_to_index[s_child]] += 0.25
    return child_distn

def get_transition_matrix(npop, sAA, sAa):
    """
    Note that sab is 0 by convention.
    @param npop: constant Wright-Fisher population
    @param sAA: a selection value
    @param sAa: a selection value
    @return: a transition matrix
    """
    fitnesses = 1.0 + np.array([sAA, sAa, 0])
    # precompute the index_to_composition and composition_to_index maps.
    compositions = list(gen_population_compositions(npop))
    c_to_i = dict((c, i) for i, c in enumerate(compositions))
    nstates = get_state_space_size(npop)
    if nstates != len(compositions):
        raise ValueError('internal error regarding state space size')
    #
    P = np.zeros((nstates, nstates))
    for parent_index, parent_composition_tuple in enumerate(compositions):
        parent_compo = np.array(parent_composition_tuple)
        random_mating = True
        if random_mating:
            single_parent_distn = parent_compo / float(np.sum(parent_compo))
            parent_distn = np.outer(single_parent_distn, single_parent_distn)
            child_distn = np.zeros(3)
            for i in range(3):
                for j in range(3):
                    child_distn += parent_distn[i, j] * get_child_distn(i, j)
            child_distn *= fitnesses
            child_distn /= np.sum(child_distn)
        else:
            total = np.dot(fitnesses, parent_compo)
            single_parent_distn = (fitnesses * parent_compo) / total
            parent_distn = np.outer(single_parent_distn, single_parent_distn)
            child_distn = np.zeros(3)
            for i in range(3):
                for j in range(3):
                    child_distn += parent_distn[i, j] * get_child_distn(i, j)
        for child_index, child_composition_tuple in enumerate(compositions):
            P[parent_index, child_index] = math.exp(
                    StatsUtil.multinomial_log_pmf(
                        child_distn, child_composition_tuple))
    return P

def get_absorbing_state_indices(npop):
    compositions = list(gen_population_compositions(npop))
    c_to_i = dict((c, i) for i, c in enumerate(compositions))
    return [
            c_to_i[(npop, 0, 0)],
            c_to_i[(0, 0, npop)],
            ]

def solve(npop, P):
    """
    @param npop: population size
    @param P: Wright-Fisher transition matrix
    @return: vector of eventual fixation probabilities of allele B
    """
    # get the absorbing state indices
    absorbing_state_indices = get_absorbing_state_indices(npop)
    # set up the system of equations
    nstates = get_state_space_size(npop)
    A = P - np.eye(nstates)
    b = np.zeros(nstates)
    for i, absorbing_state in enumerate(absorbing_state_indices):
        A[absorbing_state, absorbing_state] = 1.0
    b[absorbing_state_indices[0]] = 1.0
    # solve the system of equations
    return linalg.solve(A, b)

def G(x, c, D):
    exponent = -2*c*x*(D*(1-x) + 1)
    return math.exp(exponent)

def get_response_content(fs):
    initial_composition = (fs.nAA, fs.nAa, fs.naa)
    npop = sum(initial_composition)
    nstates = get_state_space_size(npop)
    # Check for minimum population size.
    if npop < 1:
        raise ValueError('there should be at least one individual')
    # Check the complexity;
    # solving a system of linear equations takes about n^3 effort.
    if nstates ** 3 > 1e8:
        raise ValueError('sorry this population size is too large')
    # Compute the exact probability of fixation of B.
    P = get_transition_matrix(npop, fs.sAA, fs.sAa)
    # Precompute the map from compositions to state index.
    compositions = list(gen_population_compositions(npop))
    c_to_i = dict((c, i) for i, c in enumerate(compositions))
    # Compute the exact probabilities of fixation.
    p_fixation = solve(npop, P)[c_to_i[initial_composition]]
    out = StringIO()
    print >> out, 'probability of eventual fixation (as opposed to extinction)'
    print >> out, 'of allele A in the population:'
    print >> out, p_fixation
    print >> out
    if fs.sAA:
        s = fs.sAA
        sh = fs.sAa
        h = fs.sAa / fs.sAA
        c = npop * s
        D = 2*h - 1
        p = (2*fs.nAA + fs.nAa) / float(2*npop)
        top = integrate.quad(G, 0, p, args=(c, D))[0]
        bot = integrate.quad(G, 0, 1, args=(c, D))[0]
        print >> out, 'kimura approximation:'
        print >> out, top / bot
        print >> out
    """
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
    """
    return out.getvalue()

