"""
Examine a certain rooted Laplacian matrix.

In particular, consider the Laplacian matrix
for which all internal vertices except one
have been removed by Schur complementation.
Such a matrix is related to the harmonic extension
of valuations from the leaves of a tree
to the root of the tree,
and this extension should be the same as the
Felsenstein ancestral state estimate.
"""

from StringIO import StringIO

import numpy as np
import scipy
import scipy.linalg

import Form
import FormOut
import Ftree
import FtreeIO

g_default_tree = """\
(((3:0.333333333333, 4:0.5)7:1.0, 5:1.0)8:0.088383145868,
(1:1.0, 2:0.5)6:0.911616854132);"""

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', g_default_tree)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the tree
    R, B, N = FtreeIO.newick_to_RBN(fs.tree)
    r = Ftree.R_to_root(R)
    T = Ftree.R_to_T(R)
    leaves = Ftree.T_to_leaves(T)
    internal_not_r = [v for v in Ftree.T_to_internal_vertices(T) if v is not r]
    # define the lists of leaves induced by the root
    vertex_partition = sorted(Ftree.R_to_vertex_partition(R))
    vertex_lists = [sorted(p) for p in vertex_partition]
    leaf_set = set(leaves)
    leaf_lists = [sorted(s & leaf_set) for s in vertex_partition]
    # order the list of leaves in a nice block form
    leaves = [v for lst in leaf_lists for v in lst]
    # remove internal vertices by Schur complementation
    L_schur_rooted = Ftree.TB_to_L_schur(T, B, leaves + [r])
    L_schur_full = Ftree.TB_to_L_schur(T, B, leaves)
    # show the matrix
    np.set_printoptions(linewidth=132)
    out = StringIO()
    # show the rooted schur complement
    w, v = scipy.linalg.eigh(L_schur_rooted)
    print >> out, 'rooted Schur complement:'
    print >> out, L_schur_rooted
    print >> out, 'Felsenstein weights at the root:'
    print >> out, -L_schur_rooted[-1][:-1] / L_schur_rooted[-1, -1]
    print >> out, 'rooted Schur complement eigendecomposition:'
    print >> out, w
    print >> out, v
    print >> out
    # show the full schur complement
    w, v = scipy.linalg.eigh(L_schur_full)
    print >> out, 'full Schur complement:'
    print >> out, L_schur_full
    print >> out, 'full Schur complement eigendecomposition:'
    print >> out, w
    print >> out, v
    print >> out
    # analyze perron components
    print >> out, 'perron components:'
    print >> out
    start = 0
    for lst in leaf_lists:
        n = len(lst)
        C = L_schur_rooted[start:start+n, start:start+n]
        print >> out, 'C:'
        print >> out, C
        w_eff = np.sum(C)
        b_eff = 1 / w_eff
        print >> out, 'effective conductance:'
        print >> out, w_eff
        print >> out, 'effective branch length (or resistance or variance):'
        print >> out, b_eff
        S = np.linalg.pinv(C)
        print >> out, 'C^-1 (rooted covariance-like):'
        print >> out, S
        w, v = scipy.linalg.eigh(S)
        print >> out, 'rooted covariance-like eigendecomposition:'
        print >> out, w
        print >> out, v
        print >> out, 'perron value:'
        print >> out, w[-1]
        print >> out, 'reciprocal of perron value:'
        print >> out, 1 / w[-1]
        print >> out
        start += n
    print >> out
    # analyze subtrees
    print >> out, 'subtree Laplacian analysis:'
    print >> out
    start = 0
    for lst in vertex_lists:
        n = len(lst)
        C = Ftree.TB_to_L_schur(T, B, lst + [r])
        w, v = scipy.linalg.eigh(C)
        print >> out, 'subtree Laplacian:'
        print >> out, C
        print >> out, 'eigendecomposition:'
        print >> out, w
        print >> out, v
        print >> out
        start += n
    # analyze subtrees
    print >> out, 'full Schur complement subtree analysis:'
    print >> out
    start = 0
    for lst in leaf_lists:
        n = len(lst)
        C = Ftree.TB_to_L_schur(T, B, lst + [r])
        w, v = scipy.linalg.eigh(C)
        print >> out, 'full Schur complement in subtree:'
        print >> out, C
        print >> out, 'eigendecomposition:'
        print >> out, w
        print >> out, v
        print >> out
        start += n
    return out.getvalue()

