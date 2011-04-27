"""Compute Perron values near a vertex of a tree.

Given a single specified vertex of articulation,
all other vertices of articulation in a given tree
are removed by Schur complementation.
Let L~ be the resulting Laplacian matrix relating the p+1 remaining
vertices where p is the number of leaves of the original tree.
If the specified vertex of articulation has degree d then
the vertices can be partitioned into d+1 sets
where C0 has only the specified vertex
and Ci has the original leaves, for 1<=i<=d.
The spectral norms of inverses of principal submatrices of L~
corresponding to C1, ..., Cd are the Perron values.
In other words, the Perron values are the spectral norms
of the bottleneck matrices corresponding to C1, ..., Cd.
"""

from StringIO import StringIO
from collections import defaultdict

import scipy

import Form
import FormOut
import Harmonic
import Ftree
import FtreeIO

g_default_tree = '(((3:0.333333333333, 4:0.5)7:1.0, 5:1.0)8:0.088383145868, (1:1.0, 2:0.5)6:0.911616854132)r;'

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', g_default_tree),
            Form.SingleLine('vertex',
                'distinguished vertex of articulation', 'r')]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_unique_vertex(N, name):
    """
    @param N: vertex to name
    """
    if len(set(N.values())) != len(N.values()):
        msg = 'vertex names should be unique'
        raise ValueError(msg)
    name_to_v = dict((n, v) for v, n in N.items())
    if name not in name_to_v:
        msg = 'the vertex name was not found in the newick string'
        raise ValueError(msg)
    return name_to_v[name]

def get_algebraic_connectivity(T, B, leaves):
    L_schur = Ftree.TB_to_L_schur(T, B, leaves)
    w = scipy.linalg.eigh(L_schur, eigvals_only=True)
    return w[1]

def get_response_content(fs):
    # read the tree
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    # get the distinguished vertex of articulation
    r = get_unique_vertex(N, fs.vertex)
    if r not in internal:
        msg = 'the distinguished vertex should have degree at least two'
        raise ValueError(msg)
    # Partition the leaves with respect to the given root.
    # Each set of leaves will eventually define a connected component.
    R = Ftree.T_to_R_specific(T, r)
    v_to_sinks = Ftree.R_to_v_to_sinks(R)
    # break some edges
    R_pruned = set(R)
    neighbors = Ftree.T_to_v_to_neighbors(T)[r]
    for adj in neighbors:
        R_pruned.remove((r, adj))
    T_pruned = Ftree.R_to_T(R_pruned)
    # get the leaf partition
    ordered_leaves = []
    leaf_lists = []
    for adj in neighbors:
        R_subtree = Ftree.T_to_R_specific(T_pruned, adj)
        C = sorted(b for a, b in R_subtree if b not in v_to_sinks)
        ordered_leaves.extend(C)
        leaf_lists.append(C)
    # define the vertices to keep and those to remove
    keepers = ordered_leaves + [r]
    # get the schur complement
    L_schur = Ftree.TB_to_L_schur(T, B, keepers)
    # get principal submatrices of the schur complement
    principal_matrices = []
    accum = 0
    for component_leaves in leaf_lists:
        n = len(component_leaves)
        M = L_schur[accum:accum+n, accum:accum+n]
        principal_matrices.append(M)
        accum += n
    # write the report
    out = StringIO()
    print >> out, 'algebraic connectivity:'
    print >> out, get_algebraic_connectivity(T, B, leaves)
    print >> out
    print >> out
    print >> out, 'perron values:'
    print >> out
    for M, leaf_list in zip(principal_matrices, leaf_lists):
        value = scipy.linalg.eigh(M, eigvals_only=True)[0]
        name_list = [N[v] for v in leaf_list]
        print >> out, name_list
        print >> out, value
        print >> out
    return out.getvalue()

