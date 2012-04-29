"""
Visualize steps of an exact spectral reconstruction of a tree.
"""

from StringIO import StringIO
import itertools
import math
import random

import pydot

import Form
import FormOut
import nhj
import Util
import Ftree
import FtreeIO
import NewickIO

g_default_newick = '((1:1, 2:0.5):1, (3:0.333333333333, 4:0.5):1, 5:1);'
g_felsenstein_tree_string = '(((a:0.3, b:0.1):0.25, c:0.65):0.2, (d:0.1, e:0.1):0.7);'

"""
g_upper_adjacency = np.array([
    [0, 0, 0, 0, 0, 1, 0, 0],
    [0, 0, 0, 0, 0, 2, 0, 0],
    [0, 0, 0, 0, 0, 0, 3, 0],
    [0, 0, 0, 0, 0, 0, 2, 0],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 0, 0, 0, 0]])
"""

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('newick', 'newick tree with branch lengths',
                #g_default_newick),
                #g_felsenstein_tree_string),
                NewickIO.rooted_example_tree),
            Form.CheckGroup('options', 'pre-processing options', [
                Form.CheckItem('jitter', 'jitter branch lengths')]),
            Form.RadioGroup('component_vis',
                'connected component visualization', [
                    Form.RadioItem('vis_star',
                        'an abstract star graph visualization', True),
                    Form.RadioItem('vis_complete',
                        'a less abstract complete graph visualization')]),

            ]
    return form_objects

def get_form_out():
    return FormOut.Html('reconstruction')

def get_response_content(fs):
    # read the user input
    T, B, N = FtreeIO.newick_to_TBN(fs.newick)
    # jitter branch lengths if requested
    if fs.jitter:
        mu = 0
        sigma = 1e-5
        for edge in B:
            B[edge] *= math.exp(random.gauss(0, sigma))
    # summarize the tree
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    # define the schur complement laplacian matrix
    G = Ftree.TB_to_L_schur(T, B, leaves)
    # define the fully connected schur complement graph as a Laplacian matrix
    # init the tree reconstruction state
    v_to_name = {}
    for v in leaves:
        name = N.get(v, None)
        if name is None:
            name = 'P' + chr(ord('a') + v)
        v_to_name[v] = name
    v_to_svs = dict((v, set([0])) for v in leaves)
    sv_to_vs = {0 : set(leaves)}
    edge_to_weight = {}
    for index_pair in itertools.combinations(range(nleaves), 2):
        i, j = index_pair
        leaf_pair = (leaves[i], leaves[j])
        edge_to_weight[frozenset(leaf_pair)] = -G[index_pair]
    # pairs like (-(number of vertices in supervertex sv), supervertex sv)
    active_svs = set([0])
    # initialize the sources of unique vertex and supervertex identifiers
    v_gen = itertools.count(max(leaves)+1)
    sv_gen = itertools.count(1)
    # write the output
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<body>'
    for count_pos in itertools.count(1):
        # add the graph rendering before the decomposition at this stage
        print >> out, '<div>'
        if fs.vis_star:
            print >> out, nhj.get_svg_star_components(
                    active_svs, sv_to_vs, v_to_name, v_to_svs, edge_to_weight)
        elif fs.vis_complete:
            print >> out, nhj.get_svg(
                    active_svs, sv_to_vs, v_to_name, v_to_svs, edge_to_weight)
        print >> out, '</div>'
        # update the splits
        next_active_svs = set()
        # svs can be decomposed independently in arbitrary order
        alpha_index_gen = itertools.count()
        for sv in active_svs:
            if len(sv_to_vs[sv]) > 2:
                v_new = next(v_gen)
                sv_new_a = next(sv_gen)
                sv_new_b = next(sv_gen)
                alpha_index = next(alpha_index_gen)
                alpha = chr(ord('a') + alpha_index)
                v_to_name[v_new] = 'R%s%s' % (count_pos, alpha)
                next_active_svs.add(sv_new_a)
                next_active_svs.add(sv_new_b)
                if len(sv_to_vs[sv]) == 3:
                    sv_new_c = next(sv_gen)
                    nhj.delta_wye_transform(
                            sv, v_to_svs, sv_to_vs, edge_to_weight,
                            v_new, sv_new_a, sv_new_b, sv_new_c)
                    next_active_svs.add(sv_new_c)
                else:
                    nhj.harmonic_split_transform(
                            sv, v_to_svs, sv_to_vs, edge_to_weight,
                            v_new, sv_new_a, sv_new_b)
            else:
                next_active_svs.add(sv)
        # if the set of active svs has not changed then we are done
        if active_svs == next_active_svs:
            break
        else:
            active_svs = next_active_svs
    print >> out, '</html>'
    print >> out, '</body>'
    return out.getvalue()

