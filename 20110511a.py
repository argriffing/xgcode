"""Draw an MDS with imputed internal nodes, using TikZ.
"""

import numpy as np

import Form
import FormOut
import tikz
import Ftree
import FtreeIO
import const

g_tree_string = const.read('20100730g').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                g_tree_string),
            Form.CheckGroup('check_group', 'scaling options', [
                Form.CheckItem('scale_using_eigenvalues',
                    'scale using eigenvalues', True)]),
            Form.Float('scaling_factor',
                'scaling factor', 10.0, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_vertex_line(v, x, y):
    """
    @param v: the vertex
    @param x: vertex x location
    @param y: vertex y location
    @return: tikz line
    """
    style = 'draw,shape=circle,inner sep=0pt'
    line = '\\node (%s)[%s] at (%.4f, %.4f) {};' % (v, style, x, y)
    return line

def get_edge_line(va, vb):
    """
    @param va: first vertex
    @param vb: second vertex
    @return: tikz line
    """
    line = '\\path (%s) edge node {} (%s);' % (va, vb)
    return line

def get_tikz_lines(fs):
    """
    @param fs: user input
    @return: a sequence of tikz lines
    """
    newick = fs.tree_string
    # hardcode the axes
    x_index = 0
    y_index = 1
    # get the tree with ordered vertices
    #T, B = FtreeIO.newick_to_TB(newick, int)
    T, B, N = FtreeIO.newick_to_TBN(newick)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # get the harmonic extension points
    w, v = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
    # possibly scale using the eigenvalues
    if fs.scale_using_eigenvalues:
        X_full = np.dot(v, np.diag(np.reciprocal(np.sqrt(w))))
    else:
        X_full = v
    # scale using the scaling factor
    X_full *= fs.scaling_factor
    # get the first two axes
    X = np.vstack([X_full[:,x_index], X_full[:,y_index]]).T
    # get the tikz lines
    axis_lines = [
            '% draw the axes',
            '\\node (axisleft) at (0, -5) {};',
            '\\node (axisright) at (0, 5) {};',
            '\\node (axistop) at (5, 0) {};',
            '\\node (axisbottom) at (-5, 0) {};',
            '\\path (axisleft) edge[draw,color=lightgray] node {} (axisright);',
            '\\path (axistop) edge[draw,color=lightgray] node {} (axisbottom);']
    node_lines = []
    for v, (x,y) in zip(vertices, X.tolist()):
        line = get_vertex_line(v, x, y)
        node_lines.append(line)
    edge_lines = []
    for va, vb in T:
        line = get_edge_line(va, vb)
        edge_lines.append(line)
    return axis_lines + node_lines + edge_lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    tikz_body_lines = get_tikz_lines(fs)
    tikz_body = '\n'.join(tikz_body_lines)
    return tikz.get_tikz_response([], '', tikz_body, fs.tikzformat)

