"""Root the tree at the harmonic Fiedler point.

Let the harmonic zero domain be the set of points
on the edge-extended tree
where the harmonic extension of the Fiedler vector of the
Schur complement Laplacian matrix is zero.
If the harmonic zero domain consists of a single point
then we will call it the harmonic Fiedler point.
"""

from StringIO import StringIO
from collections import defaultdict

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
    return FormOut.Newick()

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
    T.remove(u_edge)
    del B[u_edge]
    ea = frozenset((r, a))
    eb = frozenset((r, b))
    T.add(ea)
    T.add(eb)
    B[ea] = da
    B[eb] = db
    R = Ftree.T_to_R_specific(T, r)
    # return the string representation of the new tree
    return FtreeIO.RBN_to_newick(R, B, N)
