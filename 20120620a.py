"""
Endpoint conditioned sampling from a population genetic model.

The model is similar to that of Wright-Fisher
and has mutation, selection, and recombination.
Each of these effects is associated with a separate global parameter.
For the initial and final state fields,
each row corresponds to a single chromosome
and each column corresponds to a position within the chromosome.
At every position in every chromosome
is either a low fitness or a high fitness allele,
denoted by 0 or 1 respectively.
"""

from StringIO import StringIO
from itertools import combinations

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
import mrate
import MatrixUtil
from smallutil import stripped_lines
import gmpy

g_default_initial_state = """
1111
0000
1010
""".strip()

g_default_final_state = """
0000
1100
0011
""".strip()

def get_form():
    form_objects = [
            Form.Integer('ngenerations', 'total number of generations',
                10, low=3, high=100),
            Form.Float('selection_param', 'multiplicative allele fitness',
                '1.01', low_inclusive=1),
            Form.Float('mutation_param', 'mutation randomization probability',
                '0.01', low_inclusive=0, high_inclusive=1),
            Form.Float('recombination_param',
                'linkage phase randomization probability',
                '0.01', low_inclusive=0, high_inclusive=1),
            Form.MultiLine('initial_state', 'initial state',
                g_default_initial_state),
            Form.MultiLine('final_state', 'final state',
                g_default_final_state),
            ]
    return form_objects

def get_form_out():
    return FormOut.Report()

def multiline_state_to_ndarray(multiline_state):
    arr = []
    for line in stripped_lines(multiline_state.splitlines()):
        row = []
        for s in line:
            v = int(s)
            if v not in (0, 1):
                raise ValueError('invalid allele')
            row.append(v)
        arr.append(row)
    return np.array(arr)

def ndarray_to_multiline_state(K):
    return '\n'.join(''.join(str(x) for x in r) for r in K)

def ndarray_to_integer_state(K):
    """
    @param K: a rectangular array of integer ones and zeros
    @return: a python integer representing the state
    """
    nchromosomes, npositions = K.shape
    w = []
    for i in range(nchromosomes):
        for j in range(npositions):
            w.append(K[i,j])
    return sum(v<<i for i, v in enumerate(w))

def integer_state_to_ndarray(k, nchromosomes, npositions):
    """
    @param k: a python integer representing the state
    @return: a rectangular array of integer ones and zeros
    """
    arr = []
    for i in range(nchromosomes):
        row = []
        for j in range(npositions):
            v = k & 1
            k >>= 1
            row.append(v)
        arr.append(row)
    return np.array(arr)

def bitphase_to_nchanges(bitphase, npositions):
    mask = (1<<(npositions-1)) - 1
    return gmpy.popcount(mask & (bitphase ^ (bitphase >> 1)))

def get_chromosome_distn(selection, recombination, K):
    """
    Define the distribution over child chromosomes.
    Given a parental population,
    the child chromosomes are independently and identically distributed
    when only selection and recombination are considered.
    This transition is independent from the immediately
    subsequent mutation effect.
    @param selection: a fitness ratio
    @param recombination: a linkage phase randomization probability
    @param K: ndarray parental population state
    @return: array defining a conditional distribution over chromosomes
    """
    nstates = 1 << (nchromosomes * npositions)
    distn = np.zeros(1<<npositions)
    # sum over all ways to independently pick parental chromosomes
    # this is (nchromosomes)^2 things because repetition is allowed
    for parent_a in range(nstates):
        weight_a = selection**np.sum(K[parent_a])
        for parent_b in range(nstates):
            weight_b = selection**np.sum(K[parent_b])
            # sum over all recombination phases
            # this is 2^(npositions) things
            parent_pairs = zip(K[parent_a], K[parent_b])
            for bitphase in range(1<<npositions):
                # count the number of phase transitions in the bitphase
                nchanges = bitphase_to_nchanges(bitphase)
                # compute the weight corresponding to this phase transition
                weight_phase = 0
                weight_phase += recombination**nchanges
                weight_phase += (1-recombination)**(npositions-1-nchanges)
                # get the corresponding chromosome index
                w = [pair[(bitphase>>i) & 1] for i, pair in enumerate(
                    parent_pairs)]
                index = sum(v<<i for i, v in enumerate(w))
                distn[index] *= weight_a * weight_b * weight_phase
    return distn / np.sum(distn)


def get_transition_matrix(
        selection, mutation, recombination,
        nchromosomes, npositions):
    """
    This factors into the product of two transition matrices.
    The first transition matrix accounts for selection and recombination.
    The second transition matrix independently accounts for mutation.
    @param selection: a fitness ratio
    @param mutation: a state randomization probability
    @param recombination: a linkage phase randomization probability
    @param nchromosomes: number of chromosomes in the population
    @param npositions: number of positions per chromosome
    @return: a numpy array
    """
    nstates = 1 << (nchromosomes * npositions)
    # init the unnormalized selection and recombination transition matrix
    M = np.zeros((nstates, nstates))
    for source_index in nstates:
        K_source = integer_state_to_ndarray(
                source_index, nchromosomes, npositions)
        chromosome_distn = get_chromosome_distn(
                selection, recombination, K_source)
        for x in product(range(1<<npositions), repeat=nchromosomes):
                weight = 1
                for index in x:
                    weight *= chromosome_distn[x]
                sink_index = 0
                for i, index in enumerate(x):
                    sink_index <<= i*npositions
                    sink_index |= index
                M[source_index, sink_index] = weight
    M_row_sums = np.sum(M, axis=1)
    P_selection_recombination = (M.T / M_row_sums).T
    # define the mutation transition matrix
    P_mutation = np.zeros((nstates, nstates))
    for source in nstates:
        for sink in states:
            d = gmpy.hamdist(source, sink)
            P_mutation[source, sink] = (mutation**d) * (1-mutation)**(nstates-d)
    # define the state transition probability matrix
    P = np.dot(P_selection_recombination, P_mutation)
    return P


def get_response_content(fs):
    initial_state = multiline_state_to_ndarray(fs.initial_state)
    final_state = multiline_state_to_ndarray(fs.final_state)
    if initial_state.shape != final_state.shape:
        raise ValueError(
                'initial and final states do not have the same dimensions')
    nchromosomes, npositions = initial_state.shape
    if nchromosomes * npositions > 12:
        raise ValueError('at most 2^12 states are allowed per generation')
    #
    out = StringIO()
    print >> out, 'number of chromosomes:'
    print >> out, nchromosomes
    print >> out
    print >> out, 'number of positions per chromosome:'
    print >> out, npositions
    print >> out
    for K in (initial_state, final_state):
        print >> out, ndarray_to_multiline_state(
                integer_state_to_ndarray(
                    ndarray_to_integer_state(K), nchromosomes, npositions))
        print >> out
    return out.getvalue()

