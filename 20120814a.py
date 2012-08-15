"""
Compute fixation probabilities for W-F models with selection and recombination.

Do this by solving a Markov chain numerically.
W-F is Wright-Fisher.
The model is assumed to be haploid,
and we consider the co-evolution of two nucleotide sites.
The type at the first nucleotide site is either 'A' or 'a',
while the type at the second nucleotide site is either 'B' or 'b'.
Note that these are two different sites on the same haploid chromosome,
as opposed to homologous types in a diploid pair.
Each of four states (AB, Ab, aB, ab) has its own fitness,
but we will arbitrarily set the fitness of 'ab' to 1.
Fitness is 1+s where s is the selection value (0 by convention for state 'ab').
Recombination probability is between 0 and 1/2;
from this range you can probably figure out how it is interpreted.
"""

from StringIO import StringIO
import math
from math import exp
import time

import numpy as np
from scipy import linalg
from scipy import interpolate

import Form
import FormOut
import MatrixUtil
import StatsUtil
import StatsVectorized
import Util
import wfengine

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('nAB', 'number of AB alleles', 1, low=0),
            Form.Integer('nAb', 'number of Ab alleles', 5, low=0),
            Form.Integer('naB', 'number of aB alleles', 6, low=0),
            Form.Integer('nab', 'number of ab alleles', 0, low=0),
            Form.Float('sAB', 'AB selection', 0.4, low_exclusive=-1),
            Form.Float('sAb', 'Ab selection', 0.2, low_exclusive=-1),
            Form.Float('saB', 'aB selection', 0.1, low_exclusive=-1),
            Form.Float('r', 'recombination probability',
                0.1, low_inclusive=0, high_inclusive=0.5),
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
    k = 4
    return Util.choose(N + k - 1, k - 1)

def gen_population_compositions(npop):
    nstates = npop + 1
    for i in range(nstates):
        for j in range(nstates - i):
            for k in range(nstates - (i+j)):
                yield (i, j, k, npop - (i+j+k))

def get_recombination(i, j):
    """
    Aggregate variants are 0:AB, 1:Ab, 2:aB, 3:ab.
    @param i: aggregate variant index
    @param j: aggregate variant index
    @return: two aggregate variant indices
    """
    index_to_string = ('AB', 'Ab', 'aB', 'ab')
    string_to_index = dict((v, i) for i, v in enumerate(index_to_string))
    si, sj = index_to_string[i], index_to_string[j]
    s_child_0 = si[0] + sj[1]
    s_child_1 = sj[0] + si[1]
    return string_to_index[s_child_0], string_to_index[s_child_1]

def get_child_distns(npop, sAB, sAb, saB, r):
    """
    Note that sab is 0 by convention.
    @param npop: constant Wright-Fisher population
    @param sAB: a selection value
    @param sAb: a selection value
    @param saB: a selection value
    @param r: a recombination probability
    @return: a list of distributions over child aggregate allele types
    """
    fitnesses = 1.0 + np.array([sAB, sAb, saB, 0])
    # precompute the index_to_composition and composition_to_index maps.
    compositions = list(gen_population_compositions(npop))
    compos = [np.array(x) for x in compositions]
    c_to_i = dict((c, i) for i, c in enumerate(compositions))
    nstates = get_state_space_size(npop)
    if nstates != len(compositions):
        raise ValueError('internal error regarding state space size')
    #
    child_distns = []
    for parent_index, parent_composition_tuple in enumerate(compositions):
        # Given the parent composition and selections,
        # we can find the distribution over parent alleles.
        # Because each parent in a pair is chosen independently,
        # this also gives the joint distribution over parent pair variants.
        # From this we can use the recombination probability
        # to get the distribution over child allele states.
        # From this child allele state distribution,
        # we can use a multinomial distribution to get the distribution
        # over variant type compositions for the next generation.
        parent_compo = np.array(parent_composition_tuple)
        single_parent_distn = parent_compo / float(np.sum(parent_compo))
        parent_distn = np.outer(single_parent_distn, single_parent_distn)
        child_distn = np.zeros(4)
        for i in range(4):
            for j in range(4):
                p_parents = parent_distn[i, j]
                child_distn[i] += p_parents * 0.5 * (1 - r)
                child_distn[j] += p_parents * 0.5 * (1 - r)
                recomb_i, recomb_j = get_recombination(i, j)
                child_distn[recomb_i] += p_parents * 0.5 * r
                child_distn[recomb_j] += p_parents * 0.5 * r
        child_distn *= fitnesses
        child_distn /= np.sum(child_distn)
        child_distns.append(child_distn)
    """
        for child_index, child_compo in enumerate(compositions):
            #P[parent_index, child_index] = StatsVectorized.multinomial_pmf(
                    #npop, child_distn, child_compo)
            P[parent_index, child_index] = math.exp(
                    StatsUtil.multinomial_log_pmf(
                        child_distn, child_compo))
    return P
    """
    return child_distns

def get_absorbing_state_indices(npop):
    compositions = list(gen_population_compositions(npop))
    c_to_i = dict((c, i) for i, c in enumerate(compositions))
    return [
            c_to_i[(npop, 0, 0, 0)],
            c_to_i[(0, npop, 0, 0)],
            c_to_i[(0, 0, npop, 0)],
            c_to_i[(0, 0, 0, npop)],
            ]

def solve(npop, P):
    """
    @param npop: population size
    @param P: Wright-Fisher transition matrix
    @return: matrix of eventual fixation probabilities
    """
    # get the absorbing state indices
    absorbing_state_indices = get_absorbing_state_indices(npop)
    # set up the system of equations
    nstates = get_state_space_size(npop)
    A = P - np.eye(nstates)
    B = np.zeros((nstates, 4))
    for i, absorbing_state in enumerate(absorbing_state_indices):
        A[absorbing_state, absorbing_state] = 1.0
        B[absorbing_state, i] = 1.0
    # solve the system of equations
    return linalg.solve(A, B)

def get_response_content(fs):
    initial_composition = (fs.nAB, fs.nAb, fs.naB, fs.nab)
    npop = sum(initial_composition)
    nstates = get_state_space_size(npop)
    # Check for minimum population size.
    if npop < 1:
        raise ValueError('there should be at least one individual')
    # Check the complexity;
    # solving a system of linear equations takes about n^3 effort.
    if npop > 32:
        raise ValueError('sorry this population size is too large')
    # Compute the transition matrix.
    tm = time.time()
    child_distns = get_child_distns(npop, fs.sAB, fs.sAb, fs.saB, fs.r)
    log_child_distns = np.log(child_distns)
    print time.time() - tm, 'seconds to compute the log child distributions'
    tm = time.time()
    log_P = wfengine.expand_multinomials(npop, log_child_distns)
    print time.time() - tm, 'seconds to compute the log transition matrix'
    tm = time.time()
    P = np.exp(log_P)
    print time.time() - tm, 'seconds to entrywise exponentiate the matrix'
    # Precompute the map from compositions to state index.
    compositions = list(gen_population_compositions(npop))
    c_to_i = dict((c, i) for i, c in enumerate(compositions))
    # Compute the exact probabilities of fixation.
    tm = time.time()
    fixation_distribution = solve(npop, P)[c_to_i[initial_composition]]
    print time.time() - tm, 'seconds to solve the linear system'
    # Write the output
    out = StringIO()
    print >> out, 'distribution over eventual allele fixations:'
    for p, name in zip(fixation_distribution, ('AB', 'Ab', 'aB', 'ab')):
        print >> out, name, ':', p
    """
    print >> out, 'probability of eventual fixation (as opposed to extinction)'
    print >> out, 'of allele B in the population:'
    print >> out, p_fixation
    print >> out
    print >> out, 'Kimura would give the approximation'
    print >> out, k_top / k_bot
    print >> out
    """
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

