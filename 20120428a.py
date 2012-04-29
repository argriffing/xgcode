"""
Visualize steps of a spectral reconstruction of a graph approximating a tree.

Independent random exponentially distributed errors are added
to each edge weight of the Schur complement Laplacian matrix
before tree reconstruction is attempted.
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

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('newick', 'newick tree with branch lengths',
                #g_default_newick),
                #g_felsenstein_tree_string),
                NewickIO.rooted_example_tree),
            Form.Float('weight_delta_mu',
                'add random weight errors with this mean',
                '0.01', low_exclusive=0),
            ]
    return form_objects

def get_form_out():
    return FormOut.Html('reconstruction')

def get_svg(active_svs, sv_to_vs, v_to_name, v_to_svs, edge_to_weight):
    # get the set of active vertices
    vs = set()
    for sv in active_svs:
        vs.update(sv_to_vs[sv])
    # initialize the graph
    pydot_graph = pydot.Dot(
            graph_type='graph',
            overlap='0',
            sep='0.01',
            size='8, 8',
            )
    # define the pydot node objects with the right names
    v_to_pydot_node = dict((v, pydot.Node(v_to_name[v])) for v in vs)
    # define the edges
    edge_to_pydot_edge = {}
    for sv in active_svs:
        nvertices = len(sv_to_vs[sv])
        for pair in itertools.combinations(sv_to_vs[sv], 2):
            edge = frozenset(pair)
            distance = 1 / edge_to_weight[edge]
            va, vb = pair
            pna = v_to_pydot_node[va]
            pnb = v_to_pydot_node[vb]
            if nvertices == 2:
                # annotate the edge with the branch length
                label = '%.3f' % distance
                pydot_edge = pydot.Edge(pna, pnb, label=label)
            else:
                # do not annotate the edge with the branch length
                pydot_edge = pydot.Edge(pna, pnb)
            edge_to_pydot_edge[edge] = pydot_edge
    # add the nodes and edges
    for pydot_node in v_to_pydot_node.values():
        pydot_graph.add_node(pydot_node)
    for pydot_edge in edge_to_pydot_edge.values():
        pydot_graph.add_edge(pydot_edge)
    # do the physical layout and create the svg string
    tmp_path = Util.create_tmp_file(data=None, prefix='tmp', suffix='.svg')
    pydot_graph.write_svg(tmp_path, prog='neato')
    with open(tmp_path) as fin:
        svg_str = fin.read()
    # return the svg except for the first few lines
    svg_str = '\n'.join(svg_str.splitlines()[6:])
    return svg_str

def get_response_content(fs):
    # read the user input
    weight_delta_mu = fs.weight_delta_mu
    T, B, N = FtreeIO.newick_to_TBN(fs.newick)
    # summarize the tree
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    # define the schur complement laplacian matrix
    G = Ftree.TB_to_L_schur(T, B, leaves)
    # add some random amount to each edge weight
    for i in range(nleaves):
        for j in range(i):
            rate = 1 / fs.weight_delta_mu
            x = random.expovariate(rate)
            G[i, j] -= x
            G[j, i] -= x
            G[i, i] += x
            G[j, j] += x
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
        print >> out, get_svg(
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

