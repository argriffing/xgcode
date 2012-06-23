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
from itertools import combinations, product
import random

import numpy as np
from numpy import linalg
import gmpy

import Form
import FormOut
import MatrixUtil
from smallutil import stripped_lines
import Util
import popgenmarkov
import pgmfancy

g_default_initial_state = """
1111
0000
""".strip()

g_default_final_state = """
1100
0011
""".strip()

def get_form():
    form_objects = [
            Form.Integer('ngenerations', 'total number of generations',
                10, low=3, high=20),
            Form.Float('selection_param', 'multiplicative allele fitness',
                '2.0', low_exclusive=0),
            Form.Float('mutation_param', 'mutation randomization probability',
                '0.0001', low_inclusive=0, high_inclusive=1),
            Form.Float('recombination_param',
                'linkage phase randomization probability',
                '0.001', low_inclusive=0, high_inclusive=1),
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

#TODO obsolete remove
def get_transition_matrix(
        selection, mutation, recombination,
        nchromosomes, npositions):
    """
    This factors into the product of two transition matrices.
    The first transition matrix accounts for selection and recombination.
    The second transition matrix independently accounts for mutation.
    @param selection: a fitness ratio
    @param mutation: a state change probability
    @param recombination: a linkage phase change probability
    @param nchromosomes: number of chromosomes in the population
    @param npositions: number of positions per chromosome
    @return: a numpy array
    """
    P_mutation = popgenmarkov.get_mutation_transition_matrix(
            mutation, nchromosomes, npositions)
    P_selection_recombination = popgenmarkov.get_selection_recombination_transition_matrix(
            selection, recombination, nchromosomes, npositions)
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
    """
    # create transition matrices raised to various powers
    max_power = path_length - 2
    matrix_powers = [1]
    matrix_powers.append(P)
    for i in range(2, max_power + 1):
        matrix_powers.append(np.dot(matrix_powers[i-1], P))
    """
    # sample the path
    path = [initial_state]
    for i in range(1, path_length-1):
        M = linalg.matrix_power(P, path_length - i)
        previous_state = path[i-1]
        weight_state_pairs = []
        for state_index in range(nstates):
            weight = 1
            weight *= P[previous_state, state_index]
            weight *= M[state_index, final_state]
            #weight *= matrix_powers[path_length - 1 - i][state_index, final_state]
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
    ndimcap = 12
    if nchromosomes * npositions > ndimcap:
        raise ValueError(
                'at most 2^%d states are allowed per generation' % ndimcap)
    #
    mutation = 0.5 * fs.mutation_param
    recombination = 0.5 * fs.recombination_param
    selection = fs.selection_param
    #
    out = StringIO()
    print >> out, 'number of chromosomes:'
    print >> out, nchromosomes
    print >> out
    print >> out, 'number of positions per chromosome:'
    print >> out, npositions
    print >> out
    # define the transition matrix
    results = pgmfancy.get_state_space_info(nchromosomes, npositions)
    ci_to_short, short_to_count, sorted_chrom_lists = results
    initial_ci = pgmfancy.chroms_to_index(
            sorted(popgenmarkov.bin_to_int(row) for row in initial_state),
            npositions)
    initial_short = ci_to_short[initial_ci]
    final_ci = pgmfancy.chroms_to_index(
            sorted(popgenmarkov.bin_to_int(row) for row in final_state),
            npositions)
    final_short = ci_to_short[final_ci]
    P_sr_s = pgmfancy.get_selection_recombination_transition_matrix_s(
            ci_to_short, short_to_count, sorted_chrom_lists,
            selection, recombination, nchromosomes, npositions)
    P_m_s = pgmfancy.get_mutation_transition_matrix_s(
            ci_to_short, short_to_count, sorted_chrom_lists,
            mutation, nchromosomes, npositions)
    P = np.dot(P_sr_s, P_m_s)
    # sample the endpoint conditioned path
    path = sample_endpoint_conditioned_path(
            initial_short, final_short,
            fs.ngenerations, P)
    print >> out, 'sampled endpoint conditioned path, including endpoints:'
    print >> out
    # print the initial state without shuffling of individuals
    print >> out, ndarray_to_multiline_state(initial_state)
    print >> out
    # print the intermediate states with randomly permuted individuals
    for short_index in path[1:-1]:
        chroms = sorted_chrom_lists[short_index]
        random.shuffle(chroms)
        K = np.array([popgenmarkov.int_to_bin(x, npositions) for x in chroms])
        print >> out, ndarray_to_multiline_state(K)
        print >> out
    # print the final state without shuffling of individuals
    print >> out, ndarray_to_multiline_state(final_state)
    print >> out
    return out.getvalue()

