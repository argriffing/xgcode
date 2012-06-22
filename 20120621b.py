"""
Compute some endpoint-conditioned probabilities of properties of a path.

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

import numpy as np
from numpy import linalg

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
                10, low=2, high=100),
            Form.Float('selection_param', 'multiplicative allele fitness',
                '2.0', low_exclusive=0),
            Form.Float('mutation_param',
                'mutation probability per opportunity',
                #'per site per individual per generation boundary',
                '0.0001', low_inclusive=0, high_inclusive=1),
            Form.Float('recombination_param',
                'linkage phase change probability per opportunity',
                #'per site boundary per individual per generation boundary',
                '0.001', low_inclusive=0, high_inclusive=1),
            Form.MultiLine('initial_state', 'initial state',
                g_default_initial_state),
            Form.MultiLine('final_state', 'final state',
                g_default_final_state),
            ]
    return form_objects

def get_presets():
    return [
            Form.Preset(
                'a mutation probably happened',
                {
                    'ngenerations' : '10',
                    'selection_param' : '2.0',
                    'mutation_param' : '0.0001',
                    'recombination_param' : '0.001',
                    'initial_state' : '1111\n0000\n',
                    'final_state' : '0101\n0101\n'}),
            Form.Preset(
                'a mutation must have happened',
                {
                    'ngenerations' : '10',
                    'selection_param' : '2.0',
                    'mutation_param' : '0.0001',
                    'recombination_param' : '0.001',
                    'initial_state' : '1111\n0001\n',
                    'final_state' : '1100\n0011\n'}),
            Form.Preset(
                'identical endpoints',
                {
                    'ngenerations' : '10',
                    'selection_param' : '2.0',
                    'mutation_param' : '0.0001',
                    'recombination_param' : '0.001',
                    'initial_state' : '1111\n1111\n',
                    'final_state' : '1111\n1111\n'}),
                ]

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

def get_response_content(fs):
    initial_state = multiline_state_to_ndarray(fs.initial_state)
    final_state = multiline_state_to_ndarray(fs.final_state)
    if initial_state.shape != final_state.shape:
        raise ValueError(
                'initial and final states do not have the same dimensions')
    nchromosomes, npositions = initial_state.shape
    if nchromosomes < 2:
        raise ValueError(
                'please use a population of more than one chromosome')
    if npositions < 2:
        raise ValueError(
                'please use chromosomes longer than one position')
    ndimcap = 12
    if nchromosomes * npositions > ndimcap:
        raise ValueError(
                'at most 2^%d states are allowed per generation' % ndimcap)
    # define the initial and final integer states
    initial_integer_state = popgenmarkov.bin2d_to_int(initial_state)
    final_integer_state = popgenmarkov.bin2d_to_int(final_state)
    #
    mutation = fs.mutation_param
    recombination = fs.recombination_param
    selection = fs.selection_param
    # define the prior probabilities of properties of the history
    nsiteboundaries = npositions - 1
    ngenboundaries = fs.ngenerations - 1
    no_mutation_prior = (1 - mutation)**(
            npositions*ngenboundaries*nchromosomes)
    no_recombination_prior = (1 - recombination)**(
            nsiteboundaries*ngenboundaries*nchromosomes)
    # get the transition matrices
    P_sr = popgenmarkov.get_selection_recombination_transition_matrix(
            selection, recombination, nchromosomes, npositions)
    P_s = pgmfancy.get_selection_transition_matrix(
            selection, nchromosomes, npositions)
    P_m = popgenmarkov.get_mutation_transition_matrix(
            mutation, nchromosomes, npositions)
    # define some conditional probabilities
    p_b_given_a = linalg.matrix_power(
            np.dot(P_sr, P_m), fs.ngenerations-1)[
                    initial_integer_state, final_integer_state]
    p_b_given_a_no_mutation = linalg.matrix_power(
            P_sr, fs.ngenerations-1)[
                    initial_integer_state, final_integer_state]
    p_b_given_a_no_recombination = linalg.matrix_power(
            np.dot(P_s, P_m), fs.ngenerations-1)[
                    initial_integer_state, final_integer_state]
    # define the conditional properties of properties of the history
    no_mutation_posterior = (
            no_mutation_prior * p_b_given_a_no_mutation) / (
                    p_b_given_a)
    no_recombination_posterior = (
            no_recombination_prior * p_b_given_a_no_recombination) / (
                    p_b_given_a)
    #
    out = StringIO()
    print >> out, 'probability of no mutation in the path history:'
    print >> out, 'unconditional:        %s' % no_mutation_prior
    print >> out, 'endpoint-conditioned: %s' % no_mutation_posterior
    print >> out
    print >> out, 'probability of no recombination in the path history:'
    print >> out, 'unconditional:        %s' % no_recombination_prior
    print >> out, 'endpoint-conditioned: %s' % no_recombination_posterior
    print >> out
    return out.getvalue()

