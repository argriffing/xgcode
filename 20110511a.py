"""Draw an MDS with imputed internal nodes, using TikZ.
"""


from StringIO import StringIO
import math

import numpy as np

import Form
import FormOut
import MatrixUtil
import EigUtil
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
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_text(tikz_body):
    tikz_header = '\\begin{tikzpicture}[auto]'
    tikz_footer = '\\end{tikzpicture}'
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_latex_text(tikz_text):
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\begin{document}'])
    latex_body = tikz_text
    latex_footer = '\\end{document}'
    return '\n'.join([latex_header, latex_body, latex_footer])

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

def get_tikz_lines(newick):
    """
    @param newick: a newick tree string
    @return: a sequence of tikz lines
    """
    # hardcode the axes
    x_index = 0
    y_index = 1
    # get the tree with ordered vertices
    T, B = FtreeIO.newick_to_TB(newick, int)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # get the harmonic extension points
    w, v = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
    # do not scale using eigenvalues!
    #X_full = np.dot(v, np.diag(np.reciprocal(np.sqrt(w))))
    X_full = v
    X = np.vstack([X_full[:,x_index], X_full[:,y_index]]).T
    # get the tikz lines
    axis_lines = [
            '% draw the axes',
            '\\node (axisleft) at (0, -1.2) {};',
            '\\node (axisright) at (0, 1.2) {};',
            '\\node (axistop) at (1.2, 0) {};',
            '\\node (axisbottom) at (-1.2, 0) {};',
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
    # get the texts
    tikz_lines = get_tikz_lines(fs.tree_string)
    tikz_text = get_tikz_text('\n'.join(tikz_lines))
    latex_text = get_latex_text(tikz_text)
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)
