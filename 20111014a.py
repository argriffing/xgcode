"""
Draw a reduced dimension commute time matrix of a codon rate matrix.

Use a continuous time symmetric Markov model of changes
among non-stop codons.
If two non-stop codons differ by a transition
then the symmetric rate is alpha.
If two non-stop codons differ by a transversion
then the symmetric rate is beta.
If two non-stop codons differ by more than one nucleotide
then the symmetric rate is zero.
Therefore the underlying DNA model is like the K80 Kimura model.
Because the rate matrix is symmetric,
the stationary distribution is uniform over all 61 non-stop codons.
Tryptophan is marked with an olive colored dot.
In this two dimensional reduction,
some codons share the same point.
This sharing is caused by the symmetry of the genetic code,
even after having reduced the symmetry by removing the stop codons
and differentiating between transitions and transversions.
"""

from StringIO import StringIO

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from Codon import g_sorted_non_stop_codons
from Codon import g_stop_codons
import tikz
import latexutil
import color

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Float('ts_rate', 'transition rate', '3.0'),
            Form.Float('tv_rate', 'transversion rate', '2.0'),
            Form.TikzFormat()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def is_nt_transition(nta, ntb):
    transitions = (
            ('A', 'T'),
            ('T', 'A'),
            ('C', 'G'),
            ('G', 'C'))
    return (nta, ntb) in transitions

def is_nt_transversion(nta, ntb):
    transversions = (
            ('A', 'C'),
            ('A', 'G'),
            ('T', 'C'),
            ('T', 'G'),
            ('C', 'A'),
            ('C', 'T'),
            ('G', 'A'),
            ('G', 'T'))
    return (nta, ntb) in transversions

def is_nt_identity(nta, ntb):
    return nta == ntb

def get_tikz_lines(fs):
    codons = g_sorted_non_stop_codons
    codon_to_index = dict((c, i) for i, c in enumerate(codons))
    # Compute all codon transitions and all codon transversions
    # among non-stop codons, including both directions.
    all_codon_transitions = []
    all_codon_transversions = []
    for source_codon in codons:
        source_nts = list(source_codon)
        for i in range(3):
            source_nt = source_nts[i]
            for sink_nt in 'ACGT':
                sink_nts = source_nts[:]
                sink_nts[i] = sink_nt
                sink_codon = ''.join(sink_nts)
                if sink_codon not in g_stop_codons:
                    codon_pair = (source_codon, sink_codon)
                    if is_nt_transition(source_nt, sink_nt):
                        all_codon_transitions.append(codon_pair)
                    elif is_nt_transversion(source_nt, sink_nt):
                        all_codon_transversions.append(codon_pair)
    # create the weighted adjacency matrix
    A = np.zeros((61, 61))
    for source_codon, sink_codon in all_codon_transitions:
        source_index = codon_to_index[source_codon]
        sink_index = codon_to_index[sink_codon]
        A[source_index, sink_index] = fs.ts_rate
    for source_codon, sink_codon in all_codon_transversions:
        source_index = codon_to_index[source_codon]
        sink_index = codon_to_index[sink_codon]
        A[source_index, sink_index] = fs.tv_rate
    # create the Laplacian matrix
    L = np.diag(np.sum(A, axis=0)) - A
    # get a chunk of the eigendecomposition of L
    w, v = scipy.linalg.eigh(L, eigvals=(1, 2))
    # get MDS points as rows of X
    X = np.dot(v, np.diag(np.sqrt(np.reciprocal(w))))
    # init the lines
    lines = []
    # draw the edges with colors according to transition or transversion
    X *= 20
    for source_codon, sink_codon in all_codon_transitions:
        source_index = codon_to_index[source_codon]
        sink_index = codon_to_index[sink_codon]
        if source_index < sink_index:
            line = '\\draw[color=w-blue] (%0.4f,%0.4f) -- (%0.4f,%0.4f);' % (
                    X[source_index, 0], X[source_index, 1],
                    X[sink_index, 0], X[sink_index, 1])
            lines.append(line)
    for source_codon, sink_codon in all_codon_transversions:
        source_index = codon_to_index[source_codon]
        sink_index = codon_to_index[sink_codon]
        if source_index < sink_index:
            line = '\\draw[color=w-red] (%0.4f,%0.4f) -- (%0.4f,%0.4f);' % (
                    X[source_index, 0], X[source_index, 1],
                    X[sink_index, 0], X[sink_index, 1])
            lines.append(line)
    # draw a dot onto tryptophan
    tryptophan_index = codon_to_index['TGG']
    line = '\\fill[color=w-olive] (%0.4f,%0.4f) circle (0.1em);' % (
            X[tryptophan_index, 0], X[tryptophan_index, 1])
    lines.append(line)
    # return the lines
    return lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body = '\n'.join(get_tikz_lines(fs))
    tikzpicture = tikz.get_picture(tikz_body, 'auto')
    return tikz.get_response(
            tikzpicture, fs.tikzformat,
            tikz.get_w_color_package_set(), tikz.get_w_color_preamble())

