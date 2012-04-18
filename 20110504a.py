"""Make more TikZ tree figures.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import tikz
import latexutil
import SpatialTree
import Newick
import Ftree
import FtreeIO
import FastDaylightLayout
import EqualArcLayout

g_tree_a = '((1:1, 2:1)x:1, 3:1, 4:1)y;'
g_tree_b = '((1:1, 2:1)x:1, (3:1, 4:1)y:1, (5:1, 6:1)z:1)w;'

g_figure_label = 'fig:with-and-without-deep-vertices'
g_figure_caption = """\
In tree (a) every internal vertex is adjacent to a leaf,
whereas in tree (b) the deepest vertex of articulation
is not adjacent to any leaf."""

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('scaling_factor', 'scaling factor',
                1.0, low_exclusive=0),
            Form.LatexFormat()]
    return form_objects

def get_form_out():
    return FormOut.Latex()

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
    options = {
            'inner sep' : '0pt',
            'scale' : scaling_factor}
    tikz_body = '\n'.join(get_tikz_lines(newick))
    return tikz.get_picture(tikz_body, **options)

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

def get_latex_response(scaling_factor, latexformat):
    figure_body_lines = [
        '\\subfloat[]{\\label{fig:without-deep-vertex}',
        get_tikz_text(g_tree_a, scaling_factor),
        '}',
        '\\subfloat[]{\\label{fig:with-deep-vertex}',
        get_tikz_text(g_tree_b, scaling_factor),
        '}']
    figure_body = '\n'.join(figure_body_lines)
    packages = ['tikz', 'subfig']
    preamble = ''
    return latexutil.get_centered_figure_response(
            figure_body, latexformat, g_figure_caption, g_figure_label,
            packages, preamble)

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    return get_latex_response(fs.scaling_factor, fs.latexformat)

def main(options):
    print get_latex_response(1.0, fs.latexformat)

if __name__ == '__main__':
    main()

