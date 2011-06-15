"""
Draw a TikZ MDS plus clouds of MDS error points.
"""

import random

import numpy as np
import scipy

import Form
import FormOut
import tikz
import Ftree
import FtreeIO
import MatrixUtil
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
            #Form.CheckGroup('check_group', 'scaling options', [
                #Form.CheckItem('scale_using_eigenvalues',
                    #'scale using eigenvalues', True)]),
            Form.Integer('nsamples',
                'use this many perturbed distance matrices',
                10, low=0, high=100),
            Form.Float('stddev',
                'standard deviation of distance errors',
                0.1, low_exclusive=0),
            Form.Float('scaling_factor',
                'scaling factor', 5.0, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

class FigureInfo:
    """
    This class aggregates info about the figure.
    Input is a tree with branch lengths,
    the number of random matrices,
    and the entrywise variance.
    From this input we construct the MDS tips of the tree using the true tree,
    the internal vertices of the tree using the harmonic extension,
    the MDS tips of distance matrices with error,
    being careful to reflect the projections towards the original points.
    """
    def __init__(self, T, B):
        """
        @param T: tree topology
        @param B: branch lengths
        """
        self.T = T
        self.B = B
        # define a leaf and internal vertex order
        self.leaves = Ftree.T_to_leaves(self.T)
        self.internal = Ftree.T_to_internal_vertices(self.T)
        # compute the reference MDS for canonical reflection
        D = Ftree.TB_to_D(self.T, self.B, self.leaves)
        self.reference_MDS = self._D_to_MDS(D)
    def _get_D_perturbed(self, D, stddev):
        """
        @param D: distance matrix as a numpy array
        @param stddev: standard deviation of perturbations of entries
        @return: a perturbed distance matrix
        """
        D_perturbed = np.array(D)
        n = len(D)
        for i in range(n):
            for j in range(i+1, n):
                e = random.gauss(0, stddev)
                D_perturbed[i, j] += e
                D_perturbed[j, i] += e
        return D_perturbed
    def _get_D(self):
        """
        @return: distance matrix relating leaves
        """
        return Ftree.TB_to_D(self.T, self.B, self.leaves)
    def _D_to_MDS(self, D):
        """
        Let the MDS be an Nx2 matrix of scaled points.
        The points should be appropriately scaled by eigenvalues.
        The reflections (column signs) of the points are arbitrary.
        The points should be ordered conformantly with the NxN distance matrix.
        @param D: NxN distance matrix
        @return: Nx2 MDS matrix
        """
        G_neg = MatrixUtil.double_centered(D) / 2.0
        # Get eigenvalues in increasing order as W,
        # and get corresponding eigenvectors as columns of V.
        W_neg, V = scipy.linalg.eigh(G_neg, eigvals=(0,1))
        # Scale columns by square root of negative eigenvalue.
        # Use absolute value to hackily deal with non-Euclidean
        # distance matrices which may arise as a result of
        # adding a perturbation.
        return V * np.sqrt(np.abs(W_neg))
    def _reflected_to_reference(self, MDS):
        """
        Return a new matrix.
        The new matrix will have columns scaled by +/- towards the reference.
        @param MDS: Nx2 MDS matrix
        """
        # compute signs of dot products of corresponding columns
        signs = []
        for a, b in zip(MDS.T, self.reference_MDS.T):
            signs.append(1 if np.dot(a, b) > 0 else -1)
        # return a copy of the MDS with the new signs
        return np.array(MDS) * signs
    def gen_point_samples(self, nsamples, stddev):
        """
        Yield maps from a vertex to a point.
        @param nsamples: sample this many perturbed distance matrices
        @param stddev: standard deviation of perturbations of entries
        """
        D = self._get_D()
        for sample_index in range(nsamples):
            D_perturbed = self._get_D_perturbed(D, stddev)
            MDS = self._D_to_MDS(D_perturbed)
            MDS = self._reflected_to_reference(MDS)
            yield dict((self.leaves[i], tuple(pt)) for i, pt in enumerate(MDS))
    def get_v_to_point(self):
        """
        This uses the harmonic extension.
        Also it uses the reference MDS for reflection.
        @return: a map from vertex to point
        """
        Lbb = Ftree.TB_to_L_block(self.T, self.B, self.internal, self.internal)
        Lba = Ftree.TB_to_L_block(self.T, self.B, self.internal, self.leaves)
        L_schur = Ftree.TB_to_L_schur(self.T, self.B, self.leaves)
        W, V = scipy.linalg.eigh(L_schur, eigvals=(1, 2))
        V = V * np.reciprocal(np.sqrt(W))
        V = self._reflected_to_reference(V)
        Y = -np.dot(np.dot(np.linalg.pinv(Lbb), Lba), V)
        MDS = np.vstack([V, Y])
        vertices = self.leaves + self.internal
        return dict((vertices[i], tuple(pt)) for i, pt in enumerate(MDS))


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

def get_point_line(x, y):
    """
    This is like an anonymous vertex.
    @return: tikz line
    """
    style = 'draw,shape=circle,inner sep=0pt'
    line = '\\node[%s] at (%.4f, %.4f) {};' % (style, x, y)
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
    T, B, N = FtreeIO.newick_to_TBN(newick)
    # get the tree with ordered vertices
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # get the tikz lines
    axis_lines = [
            '% draw the axes',
            '\\node (axisleft) at (0, -6) {};',
            '\\node (axisright) at (0, 6) {};',
            '\\node (axistop) at (6, 0) {};',
            '\\node (axisbottom) at (-6, 0) {};',
            '\\path (axisleft) edge[draw,color=lightgray] node {} (axisright);',
            '\\path (axistop) edge[draw,color=lightgray] node {} (axisbottom);']
    # set up the figure info
    info = FigureInfo(T, B)
    nsamples = 10
    # define the points caused by MDS of distance matrices with errors
    point_lines = []
    for v_to_point in info.gen_point_samples(fs.nsamples, fs.stddev):
        for x, y in v_to_point.values():
            line = get_point_line(fs.scaling_factor*x, fs.scaling_factor*y)
            point_lines.append(line)
    # get the tikz corresponding to the tree drawn inside the MDS plot
    node_lines = []
    for v, (x,y) in info.get_v_to_point().items():
        line = get_vertex_line(v, fs.scaling_factor*x, fs.scaling_factor*y)
        node_lines.append(line)
    edge_lines = []
    for va, vb in T:
        line = get_edge_line(va, vb)
        edge_lines.append(line)
    # return the tikz
    return axis_lines + point_lines + node_lines + edge_lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_lines = get_tikz_lines(fs)
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

