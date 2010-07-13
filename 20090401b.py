"""Given some sequences, find vectors whose correlations are match proportions.

Given strings, find vectors whose correlations are
proportions of matching sites.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import MatrixUtil
import iterutils
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default sequences
    default_sequences = [
            'acgt',
            'acga',
            'acaa']
    # define the list of form objects
    form_objects = [
            Form.MultiLine('sequences', 'one sequence per line',
                '\n'.join(default_sequences)),
            Form.Float('epsilon', 'small values will be shown as zero',
                '1e-10', low_inclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_vectors(n):
    """
    The returned vectors have some interesting properties.
    The dot product of each vector with the unit vector is zero.
    The dot product of each vector with itself is one.
    The dot product of each vector with each other vector is zero.
    These vectors are found using the eigendecomposition of the centering matrix.
    @param n: get this many vectors
    @return: a list of n vectors
    """
    H = np.identity(n+1) - np.ones((n+1, n+1)) / (n+1)
    w, V_t = np.linalg.eigh(H)
    V = [v.tolist() for v in V_t.T]
    sorted_wV = list(sorted((x, v) for x, v in zip(w, V)))
    sorted_w, sorted_V = zip(*sorted_wV)
    my_vectors = sorted_V[1:]
    assert len(my_vectors) == n
    assert len(my_vectors[0]) == n+1
    return my_vectors

def eps_filter(x, epsilon):
    """
    @param x: a number
    @param epsilon: a small nonnegative number
    """
    return 0 if (abs(x) < epsilon) else x

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the sequences
    sequences = []
    for raw_string in iterutils.stripped_lines(StringIO(fs.sequences)):
        sequences.append(raw_string.strip())
    # get the alphabet
    alphabet = list(sorted(set(''.join(sequences))))
    # get the vectors that should represent the symbols.
    raw_vectors = get_vectors(len(alphabet))
    # set values smaller than user-defined epsilon to zero
    vectors = [[eps_filter(x, fs.epsilon) for x in v] for v in raw_vectors]
    # map letters to vectors
    letter_to_vector = dict(zip(alphabet, vectors))
    # get the number lists corresponding to the sequences
    number_lists = []
    for sequence in sequences:
        number_list = []
        for letter in sequence:
            number_list.extend(letter_to_vector[letter])
        number_lists.append(number_list)
    # begin the response
    out = StringIO()
    # print the correlation matrix
    print >> out, MatrixUtil.m_to_string(number_lists)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
