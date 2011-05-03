"""Check the variance of the difference from the leaf mean, using math.
"""

from collections import defaultdict
from StringIO import StringIO
import random
import math

import numpy as np

import Form
import FormOut
import Harmonic
import Ftree
import FtreeIO

g_default_tree = '((1:1, 2:0.5)6:1, (3:0.333333333333, 4:0.5)7:1, 5:1)8;'

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

def quadratic_form(K, v):
    """
    This is probably not efficient.
    @param K: a 2D numpy array
    @param v: a 1D numpy array
    """
    return np.sum(np.outer(v,v) * K)

def get_new_vertex(T):
    return max(Ftree.T_to_order(T)) + 1

def add_vertex(T_in, B_in, d_edge, r, t):
    """
    Add a vertex onto an edge.
    A new tree topology and branch length assignment is returned.
    @param T_in: topology
    @param B_in: branch lengths
    @param d_edge: a directed edge
    @param r: the new vertex
    @param t: the proportion of the distance along the directed edge
    @return: a new (T, B) pair
    """
    T = set(T_in)
    B = dict(B_in)
    a, b = d_edge
    u_edge = frozenset(d_edge)
    d = B[u_edge]
    da = t*d
    db = (1-t)*d
    T.remove(u_edge)
    del B[u_edge]
    ea = frozenset((r, a))
    eb = frozenset((r, b))
    T.add(ea)
    T.add(eb)
    B[ea] = da
    B[eb] = db
    return T, B

def get_response_content(fs):
    # read the tree
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    # root arbitrarily
    R = Ftree.T_to_R_canonical(T)
    # init some sampling parameters
    npillars = 9
    # init some helper variables
    nleaves = len(leaves)
    r = get_new_vertex(T)
    vertices = internal + [r] + leaves
    combo = np.array([0]*len(internal) + [1] + [-1.0/nleaves]*nleaves)
    # Map edge position triple to the quadratic form value.
    qform = {}
    for d_edge in R:
        a, b = d_edge
        u_edge = frozenset(d_edge)
        distance = B[u_edge]
        for i in range(npillars):
            # get the proportion of the distance along the branch
            t = (i + 1) / float(npillars + 1)
            T_new, B_new = add_vertex(T, B, d_edge, r, t)
            # create the new centered covariance matrix
            L = Ftree.TB_to_L_principal(T_new, B_new, vertices)
            S = np.linalg.pinv(L)
            qform[(a, b, t*distance)] = quadratic_form(S, combo)
            #shortcombo = np.array([1] + [-1.0/nleaves]*nleaves)
            #shortvert = [r] + leaves
            #L_schur = Ftree.TB_to_L_schur(T_new, B_new, shortvert)
            #S = np.linalg.pinv(L_schur)
            #qform[(a, b, t*distance)] = quadratic_form(S, shortcombo)
    wat = sorted((val, va, vb, d) for (va, vb, d), val in qform.items())
    # write the report
    out = StringIO()
    for val, va, vb, d in wat:
        print >> out, N[va], '--[', d, ']-->', N[vb], ':', val
        print >> out
    return out.getvalue()
