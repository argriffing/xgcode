"""
Draw a tree and its harmonically extended Fiedler cut using TikZ.
"""

import scipy
from scipy import linalg
import numpy as np

import math

from SnippetUtil import HandlingError
import Newick
import FastDaylightLayout
import Form
import FormOut
import DrawEigenLacing
import tikz
import Ftree
import FtreeIO
import FtreeAux
import layout
import const 

g_tree_string = const.read('20100730g').rstrip() 


def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree = Newick.parse(g_tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Float('width', 'width (centimeters)',
                12, low_inclusive=1, high_exclusive=1000),
            Form.Float('height', 'height (centimeters)',
                8, low_inclusive=1, high_exclusive=1000),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz('tkz')

def get_latex_text(latex_body):
    return '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\usetikzlibrary{snakes}',
        '\\begin{document}',
        latex_body,
        '\\end{document}'])

def get_figure_text(figure_body):
    return '\n'.join([
        '\\begin{figure}',
        '\\centering',
        figure_body,
        #'\\caption{}',
        '\\label{fig:interlacing}',
        '\\end{figure}'])

class TikzContext:
    def __init__(self):
        self.depth = 0
        self.lines = []
        self.finished = False
        self.add_line('\\begin{tikzpicture}')
        self.depth += 1
    def add_line(self, line):
        if self.finished:
            raise ValueError('tried to add a line to a finished tikz context')
        self.lines.append((' ' * self.depth) + line)
    def finish(self):
        if not self.finished:
            self.depth -= 1
            self.add_line('\\end{tikzpicture}')
        self.finished = True
    def get_text(self):
        if not self.finished:
            raise ValueError('unfinished')
        return '\n'.join(self.lines)
    def draw_dark_line(self, x1, y1, x2, y2):
        self.add_line(
                '\\draw (%.4f,%.4f) -- (%.4f,%.4f);' % (
                    x1, y1, x2, y2))

#FIXME this is used because the analogous Ftree function does not
#      include the constant vector
def TB_to_harmonic_valuations(T, B, reflect):
    """
    @param T: topology
    @param B: branch lengths
    @return: a number of dictionaries equal to the number of leaves
    """
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    nleaves = len(leaves)
    Lbb = Ftree.TB_to_L_block(T, B, internal, internal)
    Lba = Ftree.TB_to_L_block(T, B, internal, leaves)
    L_schur = Ftree.TB_to_L_schur(T, B, leaves)
    w, v1 = scipy.linalg.eigh(L_schur)
    v2 = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), v1)
    V = np.vstack([v1, v2])
    if reflect:
        V *= -1
    vs = []
    for col in range(nleaves):
        d = dict((v, V[row, col]) for row, v in enumerate(vertices))
        vs.append(d)
    return vs

def draw_plain_branches_ftree(T, B, context, v_to_location):
    """
    Draw plain branches of a tree.
    @param T: unrooted tree topology
    @param B: branch lengths
    @param context: tikz context
    @param v_to_location: maps a vertex to a (x, y) location
    """
    # draw the multi-edges with the correct styles
    style = None
    for u_edge in T:
        pa, pb = [v_to_location[v] for v in u_edge]
        context.draw_dark_line(pa[0], pa[1], pb[0], pb[1])

def draw_labels_ftree(T, N, context, v_to_location):
    """
    Use degree anchors for label placement.
    """
    for v in Ftree.T_to_order(T):
        if v not in N:
            continue
        label = N[v]
        # get the parameters for the label
        theta = DrawEigenLacing.get_free_angle_ftree(T, v, v_to_location)
        x, y = v_to_location[v]
        # draw the text relative to the location
        theta += math.pi
        float_degree = ((theta % (2 * math.pi)) * 360) / (2 * math.pi)
        degree = int(math.floor(float_degree))
        style = 'anchor=%s,inner sep=1pt' % degree
        context.add_line(
                '\\node[%s] at (%.4f,%.4f) {%s};' % (
                    style, x, y, label))

def draw_ticks_ftree(T_in, B_in, context, v2_in, v_to_location):
    """
    Increased barb radius.
    """
    eps = 1e-8
    T = set(T_in)
    B = dict(B_in)
    v2 = dict(v2_in)
    # break the branches according to valuation signs of v2
    FtreeAux.break_branches_by_vertex_sign(T, B, v2, eps)
    # add the new locations for internal vertices
    v_to_x = dict((v, x) for v, (x, y) in v_to_location.items())
    v_to_y = dict((v, y) for v, (x, y) in v_to_location.items())
    FtreeAux.harmonically_interpolate(T, B, v_to_x)
    FtreeAux.harmonically_interpolate(T, B, v_to_y)
    vertices = Ftree.T_to_order(T)
    v_to_location = dict((v, (v_to_x[v], v_to_y[v])) for v in vertices)
    # draw the ticks
    v_to_neighbors = Ftree.T_to_v_to_neighbors(T)
    for v, location in v_to_location.items():
        x, y = location
        neighbors = v_to_neighbors[v]
        if abs(v2[v]) > eps:
            continue
        if len(neighbors) == 2:
            barb_radius_cm = 0.1
            va, vb = neighbors
            ax, ay = v_to_location[va]
            bx, by = v_to_location[vb]
            theta = math.atan2(by - ay, bx - ax)
            barbx1 = x + barb_radius_cm * math.cos(theta + math.pi/2)
            barby1 = y + barb_radius_cm * math.sin(theta + math.pi/2)
            barbx2 = x + barb_radius_cm * math.cos(theta - math.pi/2)
            barby2 = y + barb_radius_cm * math.sin(theta - math.pi/2)
            context.draw_dark_line(barbx1, barby1, barbx2, barby2)
        elif len(neighbors) > 2:
            if max(abs(v2[n]) for n in neighbors) < eps:
                continue
            context.draw_dark_dot(x, y)

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get a properly formatted newick tree with branch lengths
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    # get the vertex valuations
    reflect = False
    all_valuations = TB_to_harmonic_valuations(T, B, reflect)
    fiedler_valuations = all_valuations[1]
    # do the layout
    v_to_location = FtreeAux.equal_daylight_layout(T, B, 3)
    # get the vertex list and the initial vertex locations
    vertices = Ftree.T_to_leaves(T) + Ftree.T_to_internal_vertices(T)
    X_in = np.array([tuple(v_to_location[v]) for v in vertices])
    # fit the tree to the physical size
    physical_size = (fs.width, fs.height)
    theta = layout.get_best_angle(X_in, physical_size)
    X = layout.rotate_2d_centroid(X_in, theta)
    sz = layout.get_axis_aligned_size(X)
    sf = layout.get_scaling_factor(sz, physical_size)
    X *= sf
    # get the map from id to location for the final tree layout
    v_to_location = dict((v, tuple(r)) for v, r in zip(vertices, X))
    # draw the image
    context = TikzContext()
    draw_plain_branches_ftree(T, B, context, v_to_location)
    draw_ticks_ftree(T, B, context, fiedler_valuations, v_to_location)
    draw_labels_ftree(T, N, context, v_to_location)
    context.finish()
    tikz_text = context.get_text()
    latex_text = get_latex_text(get_figure_text(tikz_text))
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)

