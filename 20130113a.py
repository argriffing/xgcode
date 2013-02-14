"""
Check the harmonic extensions of cuts of the Schur complement in a tree.
"""

from StringIO import StringIO
import itertools
import math
import random
import functools

import numpy as np
import scipy.linalg

import Form
import FormOut
import MatrixUtil
import StoerWagner
from MatrixUtil import ndot
import combobreaker


def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.RadioGroup('cut_options', 'cut option', [
                Form.RadioItem('min_cut', 'min cut', True),
                Form.RadioItem('fiedler_cut', 'fiedler cut')]),
            Form.RadioGroup('check_options', 'compatibility check', [
                Form.RadioItem('harmonic_extension',
                    'harmonic extension', True),
                Form.RadioItem('combinatorial_extension',
                    'combinatorial extension')]),
            ]

def get_form_out():
    return FormOut.Report()

def sample_tree(nleaves):
    """
    Sample an unweighted binary tree with a given number of leaves.
    @param nleaves: number of leaves in the sampled tree
    @return: V, E
    """
    if nleaves < 3:
        raise Exception('too few requested leaves')
    V = {0, 1}
    E = {frozenset([0, 1])}
    leaves = {0, 1}
    for i in range(nleaves-2):
        v = random.choice(list(leaves))
        leaves.remove(v)
        va = len(V)
        vb = va + 1
        for nat in (va, vb):
            V.add(nat)
            leaves.add(nat)
            E.add(frozenset([v, nat]))
    return V, E

def get_unweighted_vertex_degrees(V, E):
    degrees = dict((v, 0) for v in V)
    for a, b in E:
        degrees[a] += 1
        degrees[b] += 1
    return degrees

def gen_random_weighted_binary_trees():
    """
    Generate examples.
    Maybe one of these examples will be a counterexample to some statement.
    """
    nleaves_allowed = [3, 4, 5, 6, 7, 8, 9, 10]
    #nleaves_allowed = [4]
    while True:
        ntips = random.choice(nleaves_allowed)
        narts = ntips - 2
        nverts = ntips + narts
        V, E = sample_tree(ntips)
        if len(V) != nverts:
            raise Exception('expected a nontrivial unrooted binary tree')
        # get unweighted vertex degrees
        degrees = get_unweighted_vertex_degrees(V, E)
        # get the tips and the points of articulation
        tips = sorted(v for v, d in degrees.items() if d == 1)
        arts = sorted(v for v, d in degrees.items() if d == 3)
        if len(tips) != ntips:
            raise Exception('expected a nontrivial unrooted binary tree')
        if len(arts) != narts:
            raise Exception('expected a nontrivial unrooted binary tree')
        # relabel the vertices so that leaves are first
        vtrans = dict((v, i) for i, v in enumerate(tips + arts))
        E = set(frozenset([vtrans[a], vtrans[b]]) for a, b in E)
        # construct the weighted adjacency matrix
        A = np.zeros((nverts, nverts), dtype=float)
        for a, b in E:
            w = np.random.exponential()
            A[a, b] = w
            A[b, a] = w
        yield A

def min_cut_valuator(A):
    """
    @param A: a symmetric weighted adjacency matrix
    @return: a valuation of the indices of the adjacency matrix
    """
    MatrixUtil.assert_symmetric(A)
    MatrixUtil.assert_nonnegative(A)
    MatrixUtil.assert_hollow(A)
    nverts = A.shape[0]
    v_pos = StoerWagner.stoer_wagner_min_cut(A)
    v_neg = set(range(nverts)) - v_pos
    valuations = np.zeros(nverts, dtype=float)
    for v in v_pos:
        valuations[v] = 1
    for v in v_neg:
        valuations[v] = -1
    return valuations

def fiedler_cut_valuator(A):
    """
    @param A: a symmetric weighted adjacency matrix
    @return: a valuation of the indices of the adjacency matrix
    """
    MatrixUtil.assert_symmetric(A)
    MatrixUtil.assert_nonnegative(A)
    MatrixUtil.assert_hollow(A)
    nverts = A.shape[0]
    L = np.diag(np.sum(A, axis=1)) - A
    w, v = scipy.linalg.eigh(L)
    return v[:, 1]

def harmonic_extension(A, tip_valuations):
    """
    This extension uses the 'Laplacian matrix' and the 'accompanying matrix'.
    @param A: weighted undirected adjacency matrix with tips first
    @param valuations: valuations associated with the first part of A
    @return: extension of valuations to the remaining vertices
    """
    MatrixUtil.assert_symmetric(A)
    MatrixUtil.assert_nonnegative(A)
    MatrixUtil.assert_hollow(A)
    nverts = A.shape[0]
    ntips = tip_valuations.shape[0]
    narts = nverts - ntips
    L = np.diag(np.sum(A, axis=1)) - A
    acc = -np.dot(L[:ntips, -narts:], scipy.linalg.inv(L[-narts:, -narts:]))
    if acc.shape != (ntips, narts):
        raise Exception('unexpected shape of accompanying matrix')
    MatrixUtil.assert_positive(acc)
    return np.dot(tip_valuations, acc)

def combinatorial_extension(A, tip_valuations):
    """
    This extension uses {-1, +1} valuations to minimize zero-crossings.
    It currently uses brute force.
    A more sophisticated approach would minimize not the number of
    zero-crossings, but would instead minimize the number of
    connected components induced by the valuation sign cut.
    @param A: weighted undirected adjacency matrix with tips first
    @param valuations: valuations associated with the first part of A
    @return: extension of valuations to the remaining vertices
    """
    MatrixUtil.assert_symmetric(A)
    MatrixUtil.assert_nonnegative(A)
    MatrixUtil.assert_hollow(A)
    nverts = A.shape[0]
    ntips = tip_valuations.shape[0]
    narts = nverts - ntips
    # get some arbitrarily directed edges from the adjacency matrix
    edges = []
    for i in range(nverts):
        for j in range(i+1, nverts):
            if A[i, j] > 0:
                edges.append((i, j))
    # use brute force to get the best sign valuation of internal vertices
    best_art_vals = None
    best_ncrossings = None
    for art_vals in itertools.product((-1, 1), repeat=narts):
        vfull = np.concatenate((tip_valuations, art_vals))
        ncrossings = 0
        for i, j in edges:
            if A[i, j] * vfull[i] * vfull[j] < 0:
                ncrossings += 1
        if best_art_vals is None or ncrossings < best_ncrossings:
            best_art_vals = art_vals
            best_ncrossings = ncrossings
    return np.array(best_art_vals, dtype=float)

def check_generic_cut(valuator, extendor, A):
    """
    The input matrix is expected to have a certain block structure.
    In particular, the leaf vertices are expected to
    precede the points of articulation.
    Because the tree is expected to be an unrooted binary tree,
    the relative number of leaves and points of articulation
    is determined by the size of the adjacency matrix.
    @param A: adjacency matrix of an unrooted edge-weighted binary tree
    @return: True if a counterexample is found
    """
    MatrixUtil.assert_symmetric(A)
    MatrixUtil.assert_nonnegative(A)
    MatrixUtil.assert_hollow(A)
    nverts = A.shape[0]
    if nverts < 4:
        raise Exception('expected at least four vertices')
    if nverts % 2 != 0:
        raise Exception('expected an even number of vertices')
    ntips = nverts / 2 + 1
    narts = nverts / 2 - 1
    # get the schur complement laplacian and its associated adjacency matrix
    L = np.diag(np.sum(A, axis=1)) - A
    L_tips = L[:ntips, :ntips] - ndot(
            L[:ntips, -narts:],
            scipy.linalg.inv(L[-narts:, -narts:]),
            L[-narts:, :ntips],
            )
    A_tips = np.diag(np.diag(L_tips)) - L_tips
    tip_valuations = valuator(A_tips)
    #tip_valuations -= np.mean(tip_valuations)
    #tip_valuations /= np.linalg.norm(tip_valuations)
    art_valuations = extendor(A, tip_valuations)
    #ntip_pos = sum(1 for v in tip_valuations if v > 0)
    #ntip_neg = sum(1 for v in tip_valuations if v < 0)
    #nart_pos = sum(1 for v in art_valuations if v > 0)
    #nart_neg = sum(1 for v in art_valuations if v < 0)
    #print ((ntip_pos, ntip_neg), (nart_pos, nart_neg))
    valuations = np.concatenate((tip_valuations, art_valuations))
    ncrossings = 0
    for i in range(nverts):
        for j in range(i+1, nverts):
            if valuations[i] * valuations[j] * A[i, j] < 0:
                ncrossings += 1
    if ncrossings != 1:
        # found a counterexample!
        print ncrossings
        print A
        return True


def get_response_content(fs):

    # define the amount of time we will search
    nseconds = 5

    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()

    # define the cut strategy
    if fs.min_cut:
        valuator = min_cut_valuator
    elif fs.fiedler_cut:
        valuator = fiedler_cut_valuator
    else:
        raise Exception

    # define the extension strategy
    if fs.harmonic_extension:
        extendor = harmonic_extension
    elif fs.combinatorial_extension:
        extendor = combinatorial_extension
    else:
        raise Exception

    # look for a tree for which the harmonic extension of the cut is bad
    ret = combobreaker.run_checker(
            functools.partial(check_generic_cut, valuator, extendor),
            gen_random_weighted_binary_trees(),
            nseconds=nseconds,
            niterations=None,
            )

    print >> out, ret

    # show the result
    return out.getvalue()

