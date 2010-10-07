"""Compare -(1/2)MEDE'M' to -(1/2)HMDM'H.

These values should be the same
given a whole number tip duplication factor N,
the full distance matrix D,
an 'expansion' matrix M (depends on N),
the centering matrix H,
and a weighted centering matrix E (depends on N).
"""


from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
import FelTree
import const
import MatrixUtil

g_tree_string = const.read('20100730g').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.Integer('N', 'tip duplication factor', 10)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the tree
    tree = NewickIO.parse(fs.tree_string, FelTree.NewickTree)
    # get information about the tree topology
    internal = [id(node) for node in tree.gen_internal_nodes()]
    tips = [id(node) for node in tree.gen_tips()]
    vertices = internal + tips
    ntips = len(tips)
    ninternal = len(internal)
    nvertices = len(vertices)
    # get the ordered ids with the leaves first
    ordered_ids = vertices
    # get the full distance matrix
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    # compute the two matrices to be compared
    p = ninternal
    q = ntips
    N = fs.N
    aug_a = get_aug_a(D, p, q, N)
    aug_b = get_aug_b(D, p, q, N)
    # show the output
    out = StringIO()
    print >> out, "-(1/2)MEDE'M':"
    print >> out, aug_a
    print >> out
    print >> out, "-(1/2)HMDM'H:"
    print >> out, aug_b
    print >> out
    print >> out, 'allclose:', np.allclose(aug_a, aug_b)
    return out.getvalue()

def get_M(p, q, N):
    """
    This uses the notation of Eric Stone.
    @param p: the number of internal nodes
    @param q: the number of tips
    @param N: the tip duplication factor
    """
    mshape = (p + N*q, p+q)
    a = np.hstack([np.eye(p), np.zeros((p, q))])
    b = np.hstack([np.zeros((q, p)), np.eye(q)])
    return np.vstack([a] + [b]*N)

def get_E(p, q, N):
    """
    Get a weighted centering matrix.
    This uses the notation of Eric Stone.
    @param p: the number of internal nodes
    @param q: the number of tips
    @param N: the tip duplication factor
    """
    m = np.array([1]*p + [N]*q, dtype=float)
    return np.eye(p+q) - np.outer(np.ones(p+q), m/m.sum())

def get_aug_a(D, p, q, N):
    """
    Get the -(1/2)MEDE'M' augmented matrix.
    @param D: the full distance matrix
    @param p: the number of internal nodes
    @param q: the number of tips
    @param N: the tip duplication factor
    @return: an augmented matrix for comparison
    """
    M = get_M(p, q, N)
    E = get_E(p, q, N)
    ME = np.dot(M, E)
    return -0.5 * np.dot(ME, np.dot(D, ME.T))

def get_aug_b(D, p, q, N):
    """
    Get the -(1/2)HMDM'H augmented matrix.
    @param D: the full distance matrix
    @param p: the number of internal nodes
    @param q: the number of tips
    @param N: the tip duplication factor
    @return: an augmented matrix for comparison
    """
    M = get_M(p, q, N)
    D_aug = np.dot(M, np.dot(D, M.T))
    return -0.5 * MatrixUtil.double_centered(D_aug)
