"""Make more TikZ tree figures.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import tikz
import SpatialTree
import Newick
import Ftree
import FtreeIO
import FastDaylightLayout
import EqualArcLayout

g_tree_a = '((1:1, 2:1)x:1, 3:1, 4:1)y;'
g_tree_b = '((1:1, 2:1)x:1, (3:1, 4:1)y:1, (5:1, 6:1)z:1)w;'

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('scaling_factor', 'scaling factor',
                1.0, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_vertex_line(v, depth, x, y):
    """
    @param v: the vertex
    @param depth: the vertex depth which is zero for pendent vertices
    @param x: vertex x location
    @param y: vertex y location
    @return: tikz line
    """
    if depth == 0:
        style = 'draw,shape=circle,fill=black,minimum size=3pt'
    elif depth == 1:
        style = 'draw,shape=circle,minimum size=3pt'
    else:
        style = 'draw,shape=circle'
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

def get_tikz_text(newick, scaling_factor):
    if scaling_factor != 1:
        sf = ',scale=%s' % scaling_factor
    else:
        sf = ''
    #tikz_header = r'\begin{tikzpicture}[auto%s]' % sf
    tikz_header = r'\begin{tikzpicture}[inner sep=0pt%s]' % sf
    tikz_footer = r'\end{tikzpicture}'
    tikz_body = '\n'.join(get_tikz_lines(newick))
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_tikz_lines(newick):
    tree = Newick.parse(newick, SpatialTree.SpatialTree) 
    #layout = FastDaylightLayout.StraightBranchLayout() 
    #layout.set_iteration_count(20)
    #layout.do_layout(tree) 
    EqualArcLayout.do_layout(tree)
    tree.fit((2.0, 2.0))
    name_to_location = dict((
        x.name, tree._layout_to_display(x.location)) for x in tree.preorder())
    T, B, N = FtreeIO.newick_to_TBN(newick)
    node_lines = []
    for v, depth in Ftree.T_to_v_to_centrality(T).items():
        x, y = name_to_location[N[v]]
        line = get_vertex_line(v, depth, x, y)
        node_lines.append(line)
    edge_lines = []
    for va, vb in T:
        line = get_edge_line(va, vb)
        edge_lines.append(line)
    return node_lines + edge_lines

def get_latex_text(scaling_factor):
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\usepackage{subfig}',
        '\\begin{document}'])
    figure_lines = [
        '\\begin{figure}',
        '\\centering',
        '\\subfloat[]{\\label{fig:without-deep-vertex}',
        get_tikz_text(g_tree_a, scaling_factor),
        '}',
        '\\subfloat[]{\\label{fig:with-deep-vertex}',
        get_tikz_text(g_tree_b, scaling_factor),
        '}',
        '\\caption{',
        'In tree (a) every internal vertex is adjacent to a leaf,',
        'whereas in tree (b) the deepest vertex of articulation',
        'is not adjacent to any leaf.}',
        '\\label{fig:with-and-without-deep-vertices}',
        '\\end{figure}']
    figure_text = '\n'.join(figure_lines)
    latex_footer = r'\end{document}'
    return '\n'.join([latex_header, figure_text, latex_footer])

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_text = get_tikz_text(g_tree_b, fs.scaling_factor)
    latex_text = get_latex_text(fs.scaling_factor)
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)

def main(options):
    print get_latex_text(1.0)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    options, args = parser.parse_args()
    main(options)

