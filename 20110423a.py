"""Graft a new branch such that an eigenvalue will be repeated.

The new branch will be grafted onto an existing tree
at the harmonic Fiedler point, if it exists.
This point can be determined exactly.
The length of the branch, on the other hand, will be determined numerically.
The final tree should have an algebraic connectivity of multiplicity two.
"""

from StringIO import StringIO
from collections import defaultdict

import scipy.optimize

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
            Form.MultiLine('tree', 'tree', g_default_tree),
            Form.SingleLine('root_name', 'new root name', 'r'),
            Form.SingleLine('leaf_name', 'new leaf name', 't')]
    return form_objects

def get_form_out():
    return FormOut.Newick()

def get_gap(blen, T, B, u_edge):
    """
    This function will be minimized.
    @param blen: proposed branch length
    @param T: topology
    @param B: branch lengths
    @param u_edge: undirected edge
    """
    if blen <= 0:
        return 1
    leaves = Ftree.T_to_leaves(T)
    B[u_edge] = blen
    L_schur = Ftree.TB_to_L_schur(T, B, leaves)
    ev1, ev2 = scipy.linalg.eigh(L_schur, eigvals_only=True, eigvals=(1,2))
    gap = abs(ev1 - ev2)
    return gap

def get_response_content(fs):
    # read the tree
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    # get the valuations with harmonic extensions
    w, V = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
    # get the Fiedler valuations with harmonic extensions
    h = V[:,0]
    # check for vertices with small valuations
    eps = 1e-8
    if any(abs(x)<x for x in h):
        raise ValueError('the tree has no clear harmonic Fiedler point')
    # find the edge contining the harmonic Fiedler point
    v_to_val = dict((v, h[i]) for i, v in enumerate(leaves + internal))
    d_edges = [(a,b) for a, b in T if v_to_val[a]*v_to_val[b] < 0]
    if len(d_edges) != 1:
        raise ValueError('expected the point to fall clearly on a single edge')
    d_edge = d_edges[0]
    a, b = d_edge
    # find the proportion along the directed edge
    t = v_to_val[a] / (v_to_val[a] - v_to_val[b])
    # find the distance from the new root to each endpoint vertices
    u_edge = frozenset(d_edge)
    d = B[u_edge]
    da = t*d
    db = (1-t)*d
    # create the new tree
    r = max(Ftree.T_to_order(T)) + 1
    N[r] = fs.root_name
    T.remove(u_edge)
    del B[u_edge]
    ea = frozenset((r, a))
    eb = frozenset((r, b))
    T.add(ea)
    T.add(eb)
    B[ea] = da
    B[eb] = db
    # add a new leaf with arbitrary branch length
    leaf = r + 1
    N[leaf] = fs.leaf_name
    u_edge = frozenset((r, leaf))
    T.add(u_edge)
    B[u_edge] = 1.0
    # get the best branch length to cause eigenvalue multiplicity
    blen = scipy.optimize.golden(
            get_gap, (T, B, u_edge), full_output=False, tol=1e-12)
    B[u_edge] = blen
    # return the string representation of the new tree
    R = Ftree.T_to_R_specific(T, r)
    return FtreeIO.RBN_to_newick(R, B, N)

