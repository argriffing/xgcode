"""
Sample a child chromosome given a parent population.

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
from itertools import combinations, product

import numpy as np
import gmpy

import Form
import FormOut
import MatrixUtil
from smallutil import stripped_lines
import Util

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
    nchromosomes, npositions = K.shape
    distn = np.zeros(1<<npositions)
    # sum over all ways to independently pick parental chromosomes
    # this is (nchromosomes)^2 things because repetition is allowed
    for parent_a in range(nchromosomes):
        weight_a = selection**np.sum(K[parent_a])
        for parent_b in range(nchromosomes):
            weight_b = selection**np.sum(K[parent_b])
            # sum over all recombination phases
            # this is 2^(npositions) things
            parent_pairs = zip(K[parent_a], K[parent_b])
            for bitphase in range(1<<npositions):
                # count the number of phase transitions in the bitphase
                nchanges = bitphase_to_nchanges(bitphase, npositions)
                # compute the weight corresponding to this phase transition
                weight_phase = 0
                weight_phase += recombination**nchanges
                weight_phase += (1-recombination)**(npositions-1-nchanges)
                # get the corresponding chromosome index
                w = [pair[(bitphase>>i) & 1] for i, pair in enumerate(
                    parent_pairs)]
                index = sum(v<<i for i, v in enumerate(w))
                distn[index] += weight_a * weight_b * weight_phase
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
    # define the mutation transition matrix
    P_mutation = np.zeros((nstates, nstates))
    for source in range(nstates):
        for sink in range(nstates):
            ndiff = gmpy.hamdist(source, sink)
            nsame = nchromosomes * npositions - ndiff
            P_mutation[source, sink] = (mutation**ndiff)*((1-mutation)**nsame)
    print 'mutation transition matrix row sums:'
    for value in np.sum(P_mutation, axis=1):
        print value
    print
    # init the unnormalized selection and recombination transition matrix
    M = np.zeros((nstates, nstates))
    for source_index in range(nstates):
        K_source = integer_state_to_ndarray(
                source_index, nchromosomes, npositions)
        chromosome_distn = get_chromosome_distn(
                selection, recombination, K_source)
        for x in product(range(1<<npositions), repeat=nchromosomes):
            weight = 1
            for index in x:
                weight *= chromosome_distn[index]
            sink_index = 0
            for i, index in enumerate(x):
                sink_index <<= npositions
                sink_index |= index
            M[source_index, sink_index] = weight
    M_row_sums = np.sum(M, axis=1)
    P_selection_recombination = (M.T / M_row_sums).T
    print 'selection-recombination matrix row sums:'
    for value in np.sum(P_selection_recombination, axis=1):
        print value
    print
    # define the state transition probability matrix
    P = np.dot(P_selection_recombination, P_mutation)
    return P

def sample_endpoint_conditioned_path(
        initial_state, final_state, path_length, P):
    """
    Return a sequence of states.
    The returned sequence starts at the initial state
    and ends at the final state.
    Consecutive states may be the same
    if the transition matrix has positive diagonal elements.
    This function is derived from an earlier function
    in an old SamplePath module.
    @param initial_state: the first state as a python integer
    @param final_state: the last state as a python integer
    @param path_length: the number of states in the returned path
    @param P: ndarray state transition matrix
    @return: list of integer states
    """
    # get the size of the state space and do some input validation
    MatrixUtil.assert_transition_matrix(P)
    nstates = len(P)
    if not (0 <= initial_state < nstates):
        raise ValueError('invalid initial state')
    if not (0 <= final_state < nstates):
        raise ValueError('invalid final state')
    # take care of edge cases
    if path_length == 0:
        return []
    elif path_length == 1:
        if initial_state != final_state:
            raise ValueError('unequal states for a path of length one')
        return [initial_state]
    elif path_length == 2:
        return [initial_state, final_state]
    # create transition matrices raised to various powers
    max_power = path_length - 2
    matrix_powers = [1]
    matrix_powers.append(P)
    for i in range(2, max_power + 1):
        matrix_powers.append(np.dot(matrix_powers[i-1], P))
    # sample the path
    path = [initial_state]
    for i in range(1, path_length-1):
        previous_state = path[i-1]
        weight_state_pairs = []
        for state_index in range(nstates):
            weight = 1
            weight *= P[previous_state, state_index]
            weight *= matrix_powers[path_length - 1 - i][state_index, final_state]
            weight_state_pairs.append((weight, state_index))
        next_state = Util.weighted_choice(weight_state_pairs)
        path.append(next_state)
    path.append(final_state)
    return path

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
    print >> out, 'initial and final states:'
    print >> out
    for K in (initial_state, final_state):
        print >> out, ndarray_to_multiline_state(
                integer_state_to_ndarray(
                    ndarray_to_integer_state(K), nchromosomes, npositions))
        print >> out
    print >> out
    # define the transition matrix
    P = get_transition_matrix(
            fs.selection_param, fs.mutation_param, fs.recombination_param,
            nchromosomes, npositions)
    # sample the endpoint conditioned path
    initial_integer_state = ndarray_to_integer_state(initial_state)
    final_integer_state = ndarray_to_integer_state(final_state)
    path = sample_endpoint_conditioned_path(
            initial_integer_state, final_integer_state,
            fs.ngenerations, P)
    print >> out, 'sampled endpoint conditioned path, including endpoints:'
    print >> out
    for integer_state in path:
        print >> out, ndarray_to_multiline_state(
                integer_state_to_ndarray(
                    integer_state, nchromosomes, npositions))
        print >> out
    return out.getvalue()

