"""Look at more eigendecompositions of centered distance matrices.

Use notation as in the proof we are attempting.
The block structure of the augmented distance matrix
has the internal nodes first, then the tips, then the duplicate tips.
The duplicate tips are structured by repeating blocks.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import NewickIO
import Euclid
import FelTree
import SchurAlgebra
import EigUtil
import const
import RUtil

g_tree_string = '((1:1, 2:0.5)6:1, 5:1, (3:0.2, 4:0.5)7:1)8;'

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.Integer('ndups', 'add this many extra duplicates per tip', 3),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('show_aug',
                    'show the augmented distance matrix')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_augmented_distance(D, ntips, ndups):
    """
    The full distance matrix should have a specific order.
    In particular, internal nodes should appear before tip nodes.
    @param D: the full distance matrix
    @param ntips: the number of tips of the original tree
    @param ndups: add this many extra duplicates per tip
    @return: a new distance matrix with tip extensions at the end
    """
    n = len(D)
    nnew = n + ntips*ndups
    ninternal = n - ntips
    P = np.zeros((nnew, nnew))
    for i in range(nnew):
        for j in range(nnew):
            if i < n and j < n:
                P[i,j] = D[i,j]
            elif i < n and j >= n:
                P[i,j] = D[i, ninternal + (j-n)%ntips]
            elif j < n and i >= n:
                P[i,j] = D[ninternal + (i-n)%ntips, j]
            else:
                P[i,j] = D[ninternal + (i-n)%ntips, ninternal + (j-n)%ntips]
    return P

def get_ugly_matrix(M, ninternal, ntips):
    """
    @param M: the matrix
    @param ninternal: the number of internal nodes
    @param ntips: the number of tips
    @return: a multiline string
    """
    multiline_raw = str(M)
    lines = [line.rstrip() for line in multiline_raw.splitlines()]
    lines = [x for x in lines if x]
    n = len(lines)
    if n <= ninternal + 3*ntips:
        return str(M)
    filler = ' ...,'
    out_lines = lines[:ninternal + ntips] + [filler] + lines[-2*ntips:]
    return '\n'.join(out_lines)

def get_response_content(fs):
    # build the newick tree from the string
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    ninternal = nvertices - nleaves
    # get ordered ids with the internal nodes first
    ordered_ids = get_ordered_ids(tree)
    leaf_ids = [id(node) for node in tree.gen_tips()]
    # get the distance matrix and the augmented distance matrix
    D_leaf = np.array(tree.get_partial_distance_matrix(leaf_ids))
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    D_aug = get_augmented_distance(D, nleaves, fs.ndups)
    # analyze the leaf distance matrix
    X_leaf = Euclid.edm_to_points(D_leaf)
    # get the eigendecomposition of the centered augmented distance matrix
    X_aug = Euclid.edm_to_points(D_aug, nvertices-1)
    # explicitly compute the points for the given number of dups using weights
    m = [1]*ninternal + [1+fs.ndups]*nleaves
    m = np.array(m, dtype=float) / sum(m)
    X_weighted = Euclid.edm_to_weighted_points(D, m)
    # explicitly compute the points for 10x dups
    m = [1]*ninternal + [1+fs.ndups*10]*nleaves
    m = np.array(m, dtype=float) / sum(m)
    X_weighted_10x = Euclid.edm_to_weighted_points(D, m)
    # explicitly compute the limiting points as the number of dups increases
    X = Euclid.edm_to_points(D)
    X -= np.mean(X[-nleaves:], axis=0)
    XL = X[-nleaves:]
    U, s, Vt = np.linalg.svd(XL)
    Z = np.dot(X, Vt.T)
    # report the results
    np.set_printoptions(linewidth=300, threshold=10000)
    out = StringIO()
    print >> out, 'leaf distance matrix:'
    print >> out, D_leaf
    print >> out
    print >> out, 'points derived from the leaf distance matrix'
    print >> out, '(the first column is proportional to the Fiedler vector):'
    print >> out, X_leaf
    print >> out
    if fs.show_aug:
        print >> out, 'augmented distance matrix:'
        print >> out, D_aug
        print >> out
    print >> out, 'points derived from the augmented distance matrix'
    print >> out, '(the first column is proportional to the Fiedler vector):'
    print >> out, get_ugly_matrix(X_aug, ninternal, nleaves)
    print >> out
    print >> out, 'points computed using masses:'
    print >> out, X_weighted
    print >> out
    print >> out, 'points computed using masses with 10x dups:'
    print >> out, X_weighted_10x
    print >> out
    print >> out, 'limiting points:'
    print >> out, Z
    print >> out
    return out.getvalue()

def get_ordered_ids(tree):
    """
    Maybe I could use postorder here instead.
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    return ordered_ids
