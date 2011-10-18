"""
Jointly draw some figures related to tree MDS.

Four main figures should be plotted.
First the newick tree, graphically laid out.
Second the MDS of the full distance matrix including internal vertices
and lines corresponding to edges.
Third the MDS of the leaf distance matrix
without the lines.
Fourth the MDS of the leaf distance matrix
plus harmonic extensions
plus lines corresponding to edges.
There are three label drawing modes.
The first mode is with vertex labels suppressed.
The second mode draws only the original taxon names at the leaves.
The third mode draws short vertex labels directly onto all vertices.
"""

import itertools
import string

import numpy as np
import scipy

import Form
import FormOut
import tikz
import latexutil
import Ftree
import FtreeIO
import FtreeAux
import MatrixUtil
import Newick
import const

g_tree_string = const.read('20110616a').rstrip()


def get_form():
    """
    @return: a list of form objects
    """
    # reformat the tree string
    tree = Newick.parse(g_tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 75)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree_string', 'newick tree',
                formatted_tree_string),
            Form.Float('scaling_factor',
                'figure scaling factor', 3.0, low_exclusive=0),
            Form.RadioGroup('label_mode', 'label mode', [
                Form.RadioItem('label_mode_suppressed',
                    'suppress all vertex labels', True),
                Form.RadioItem('label_mode_leaf_only',
                    'draw leaf names'),
                Form.RadioItem('label_mode_show',
                    'draw arbitrary short names at all vertices')]),
            Form.RadioGroup('sanitization_options',
                'taxon label sanitization options', [
                    Form.RadioItem('no_sanitization',
                        'allow TikZ to interpret the taxon labels', True),
                    Form.RadioItem('sanitization',
                        'attempt an ad-hoc sanitization of taxon labels')]),
            Form.RadioGroup('tree_layout', 'tree layout options', [
                Form.RadioItem('equal_arc', 'equal arc layout', True),
                Form.RadioItem('equal_daylight', 'equal daylight layout')]),
            Form.LatexFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Latex()

def gen_short_names():
    """
    Generate an unbounded number of shortish names.
    """
    for c in string.lowercase:
        yield c
    for c in string.uppercase:
        yield c
    for i in itertools.count(1):
        yield str(i)

class FigureInfo:
    """
    This class aggregates info about the figure.
    Input is a tree with branch lengths and vertex labels.
    When multiple MDS plots are involved,
    be sure to reflect the points towards the original points.
    """

    def __init__(self, T, B, N, label_mode):
        """
        @param T: tree topology
        @param B: branch lengths
        @param N: vertex to name
        @param label_mode: this is an option for MDS vertex labels
        """
        # store tree characteristics
        self.T = T
        self.B = B
        self.N = N
        # store options
        self.label_mode = label_mode
        # define a leaf and internal vertex order
        self.leaves = Ftree.T_to_leaves(self.T)
        self.internal = Ftree.T_to_internal_vertices(self.T)
        self.vertices = self.leaves + self.internal
        # map vertices to short names
        self.N_short = dict(zip(self.vertices, gen_short_names()))
        # compute the reference MDS for canonical reflection
        D_leaves = Ftree.TB_to_D(self.T, self.B, self.leaves)
        self.reference_leaf_MDS = self.get_reference_leaf_MDS()

    def get_reference_leaf_MDS(self):
        """
        Let the MDS be an Nx2 matrix of scaled points.
        The points should be appropriately scaled by eigenvalues.
        The reflections (column signs) of the points are arbitrary.
        The points should be ordered conformantly with the NxN distance matrix.
        @param D: NxN distance matrix
        @return: Nx2 MDS matrix
        """
        D_leaves = Ftree.TB_to_D(self.T, self.B, self.leaves)
        G_neg = MatrixUtil.double_centered(D_leaves) / 2.0
        # Get eigenvalues in increasing order as W,
        # and get corresponding eigenvectors as columns of V.
        W_neg, V = scipy.linalg.eigh(G_neg, eigvals=(0,1))
        # Scale columns by square root of negative eigenvalue.
        # Use absolute value to hackily deal with non-Euclidean
        # distance matrices which may arise as a result of
        # adding a perturbation.
        return V * np.sqrt(np.abs(W_neg))

    def get_full_distance_MDS(self):
        """
        Reflect towards the reference MDS.
        This is similar to the leaf MDS but for the full distance matrix.
        @return: Nx2 MDS matrix
        """
        D_full = Ftree.TB_to_D(self.T, self.B, self.vertices)
        G_neg = MatrixUtil.double_centered(D_full) / 2.0
        W_neg, V = scipy.linalg.eigh(G_neg, eigvals=(0,1))
        MDS = V * np.sqrt(np.abs(W_neg))
        # compute signs of dot products of corresponding columns
        signs = []
        nleaves = len(self.leaves)
        for a, b in zip(MDS[:nleaves].T, self.reference_leaf_MDS.T):
            signs.append(1 if np.dot(a, b) > 0 else -1)
        # return a copy of the MDS with the new signs
        return np.array(MDS) * signs

    def get_harmonically_extended_MDS(self):
        Lbb = Ftree.TB_to_L_block(self.T, self.B, self.internal, self.internal)
        Lba = Ftree.TB_to_L_block(self.T, self.B, self.internal, self.leaves)
        L_schur = Ftree.TB_to_L_schur(self.T, self.B, self.leaves)
        W, V = scipy.linalg.eigh(L_schur, eigvals=(1, 2))
        V = V * np.reciprocal(np.sqrt(W))
        V = self._reflected_to_reference(V)
        Y = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), V)
        return np.vstack([V, Y])

    def _reflected_to_reference(self, MDS):
        """
        Return a new matrix.
        The new matrix will have columns scaled by +/- towards the reference.
        @param MDS: Nx2 MDS matrix
        """
        # compute signs of dot products of corresponding columns
        signs = []
        for a, b in zip(MDS.T, self.reference_leaf_MDS.T):
            signs.append(1 if np.dot(a, b) > 0 else -1)
        # return a copy of the MDS with the new signs
        return np.array(MDS) * signs

    def _get_axis_lines(self):
        axis_radius = 0.75
        return [
                '% draw the axes',
                '\\node (aleft) at (0, %f) {};' % -axis_radius,
                '\\node (aright) at (0, %f) {};' % axis_radius,
                '\\node (atop) at (%f, 0) {};' % axis_radius,
                '\\node (abottom) at (%f, 0) {};' % -axis_radius,
                '\\path (aleft) edge[draw,color=lightgray] node {} (aright);',
                '\\path (atop) edge[draw,color=lightgray] node {} (abottom);']

    def get_tikz_tree(self, layout_mode):
        """
        Return the tikz text for the tree figure.
        @param layout_mode: equal_arc vs equal_daylight
        @return: tikz text
        """
        if layout_mode == 'equal_daylight':
            v_to_point = FtreeAux.equal_daylight_layout(self.T, self.B, 3)
        else:
            v_to_point = FtreeAux.equal_arc_layout(self.T, self.B)
        # define vertices
        vertex_lines = []
        for v, (x, y) in v_to_point.items():
            line = get_vertex_line(v, x, y)
            vertex_lines.append(line)
        # draw edges
        edge_lines = []
        for va, vb in self.T:
            line = get_edge_line(va, vb)
            edge_lines.append(line)
        # draw the labels
        label_lines = []
        if self.label_mode == 'label_mode_show':
            for v, (x, y) in v_to_point.items():
                label_lines.append(get_label_line(x, y, self.N_short[v]))
        elif self.label_mode == 'label_mode_leaf_only':
            for v in self.leaves:
                x, y = v_to_point[v]
                if v in self.N:
                    label_lines.append(get_label_line(x, y, self.N[v]))
        # return the tikz
        tikz_lines = vertex_lines + edge_lines + label_lines
        return '\n'.join(tikz_lines)

    def _get_tikz_MDS_large(self, MDS):
        """
        This can be used for any full-size MDS.
        This includes both full distance matrix MDS
        and harmonically extended leaf distance MDS.
        """
        # define axis lines
        axis_lines = self._get_axis_lines()
        # define nodes of the harmonically extended tree
        node_lines = []
        for i, v in enumerate(self.vertices):
            x, y = MDS[i]
            line = get_vertex_line(v, x, y)
            node_lines.append(line)
        # draw edges of the harmonically extended tree
        edge_lines = []
        for va, vb in self.T:
            line = get_edge_line(va, vb)
            edge_lines.append(line)
        # draw the labels
        label_lines = []
        if self.label_mode == 'label_mode_show':
            for i, v in enumerate(self.vertices):
                x, y = MDS[i]
                label_lines.append(get_label_line(x, y, self.N_short[v]))
        elif self.label_mode == 'label_mode_leaf_only':
            for i, v in enumerate(self.vertices):
                x, y = MDS[i]
                if v in self.N:
                    label_lines.append(get_label_line(x, y, self.N[v]))
        # return the tikz
        tikz_lines = axis_lines + node_lines + edge_lines + label_lines
        return '\n'.join(tikz_lines)

    def get_tikz_MDS_full(self):
        """
        Return the tikz text for the full MDS figure.
        This figure includes edges.
        @return: tikz text
        """
        MDS = self.get_full_distance_MDS()
        return self._get_tikz_MDS_large(MDS)

    def get_tikz_MDS_harmonic(self):
        """
        Return the tikz text for the harmonically extended MDS figure.
        This figure includes edges.
        @return: tikz text
        """
        MDS = self.get_harmonically_extended_MDS()
        return self._get_tikz_MDS_large(MDS)

    def get_tikz_MDS_partial(self):
        """
        Return the tikz text for the partial MDS figure.
        This figure does not include edges.
        @return: tikz text
        """
        MDS = self.reference_leaf_MDS
        # define axis lines
        axis_lines = self._get_axis_lines()
        # define leaf MDS points
        point_lines = []
        for i, v in enumerate(self.leaves):
            x, y = MDS[i]
            line = get_point_line(x, y)
            point_lines.append(line)
        # draw the labels
        label_lines = []
        if self.label_mode == 'label_mode_show':
            for i, v in enumerate(self.leaves):
                x, y = MDS[i]
                label_lines.append(get_label_line(x, y, self.N_short[v]))
        # return the tikz
        tikz_lines = axis_lines + point_lines + label_lines
        return '\n'.join(tikz_lines)


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

def get_point_line(x, y):
    """
    This is like an anonymous vertex.
    @return: tikz line
    """
    style = 'draw,shape=circle,inner sep=0pt'
    line = '\\node[%s] at (%.4f, %.4f) {};' % (style, x, y)
    return line

def get_label_line(x, y, label):
    """
    This is like an anonymous vertex.
    @return: tikz line
    """
    #style = 'draw,shape=circle,inner sep=0pt'
    style = ''
    line = '\\node[%s] at (%.4f, %.4f) {%s};' % (style, x, y, label)
    return line

def get_edge_line(va, vb):
    """
    @param va: first vertex
    @param vb: second vertex
    @return: tikz line
    """
    line = '\\path (%s) edge node {} (%s);' % (va, vb)
    return line

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    T, B, N = FtreeIO.newick_to_TBN(fs.tree_string)
    # sanitize taxon labels if requested
    if fs.sanitization:
        for v in N:
            N[v] = latexutil.sanitize(N[v])
    # scale branch lengths so the diameter is 1
    diameter = np.max(Ftree.TB_to_D(T, B, Ftree.T_to_leaves(T)))
    # scale the branch lengths
    for u_edge in T:
        B[u_edge] /= diameter
    info = FigureInfo(T, B, N, fs.label_mode)
    # get the texts
    tikz_bodies = [
            info.get_tikz_tree(fs.tree_layout),
            info.get_tikz_MDS_full(),
            info.get_tikz_MDS_partial(),
            info.get_tikz_MDS_harmonic(),
            ]
    tikz_pictures = []
    for b in tikz_bodies:
        tikzpicture = tikz.get_picture(b, 'auto', scale=fs.scaling_factor)
        tikz_pictures.append(tikzpicture)
    figure_body = '\n'.join([
        '\\subfloat[]{',
        tikz_pictures[0],
        '}',
        '\\subfloat[]{',
        tikz_pictures[1],
        '} \\\\',
        '\\subfloat[]{',
        tikz_pictures[2],
        '}',
        '\\subfloat[]{',
        tikz_pictures[3],
        '}',
        ])
    packages = ['tikz', 'subfig']
    preamble = ''
    figure_caption = None
    figure_label = None
    return latexutil.get_centered_figure_response(
            figure_body, fs.latexformat, figure_caption, figure_label,
            packages, preamble)

