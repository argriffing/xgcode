"""Make a TikZ tree figure.

This is probably useful only as an example.
"""

# Do it like this.
# python 20090617a.py > x.tex; latex x.tex; dvips x; ps2pdf x.ps; evince x.pdf

from StringIO import StringIO
import optparse
import tempfile
import subprocess
import os

from SnippetUtil import HandlingError
from Form import RadioItem
import Form


class Layout:
    """
    This is an abstract base class.
    """

    def get_tikz_text(self, scaling_factor, show_full_tree, show_pruned_trees):
        """
        @param scaling_factor: a float
        @param show_full_tree: True if the full tree should be shown
        @param show_pruned_trees: True if the pruned trees should be shown
        @return: a multiline string that defines a tikzpicture
        """
        out = StringIO()
        if scaling_factor == 1.0:
            print >> out, '\\begin{tikzpicture}[auto]'
        else:
            print >> out, '\\begin{tikzpicture}[auto,scale=%s]' % str(scaling_factor)
        print >> out, self.get_tikz_contents(show_full_tree, show_pruned_trees)
        print >> out, '\\end{tikzpicture}'
        return out.getvalue().strip()

    def get_latex_text(self, scaling_factor, show_full_tree, show_pruned_trees):
        """
        @param scaling_factor: a float
        @param show_full_tree: True if the full tree should be shown
        @param show_pruned_trees: True if the pruned trees should be shown
        @return: a multiline string that is the contents of a valid LaTeX file
        """
        out = StringIO()
        print >> out, '\\documentclass{article}'
        print >> out, '\\usepackage{tikz}'
        print >> out, '\\begin{document}'
        print >> out, self.get_tikz_text(scaling_factor, show_full_tree, show_pruned_trees)
        print >> out, '\\end{document}'
        return out.getvalue().strip()


class SixLeafLayout(Layout):

    def get_left_stuff(self):
        """
        Use a manual layout.
        @return: (left_nodes, left_edges, left_subtree_root_node)
        """
        g = Node(7, -3, 1)
        h = Node(8, -1, 0)
        a = Node(1, -3, -1)
        b = Node(2, -5, 0)
        c = Node(3, -5, 2)
        left_nodes = [a, b, c, h, g]
        left_edges = [
                (g, 4, h),
                (h, 2, a),
                (g, 2, b),
                (c, 9, g)]
        return left_nodes, left_edges, h

    def get_right_stuff(self):
        """
        Use a manual layout.
        @return: (right_nodes, right_edges, right_subtree_root_node)
        """
        i = Node(9, 3, 1)
        j = Node(10, 1, 0)
        f = Node(6, 3, -1)
        e = Node(5, 5, 0)
        d = Node(4, 5, 2)
        right_nodes = [d, e, f, i, j]
        right_edges = [
                (j, 7, i),
                (i, 1, d),
                (e, 3, i),
                (f, 2, j)]
        return right_nodes, right_edges, j

    def get_full_tree_text(self):
        """
        @return: a multiline string of TikZ code
        """
        # get the left and right subtrees
        left_nodes, left_edges, left_root = self.get_left_stuff()
        right_nodes, right_edges, right_root = self.get_right_stuff()
        # combine the subtrees
        nodes = left_nodes + right_nodes
        edges = left_edges + right_edges + [(left_root, 1, right_root)]
        return get_tree_text(nodes, edges)

    def get_left_subtree_text(self, offset):
        """
        @param offset: vertical offset
        @return: a multiline string of TikZ code
        """
        # get the left subtree
        left_nodes, left_edges, left_root = self.get_left_stuff()
        # define the pruned subtree
        right_node = Node('B', 2, 0)
        nodes = left_nodes + [right_node]
        # define the new edge
        edges = left_edges + [(left_root, '$\\frac{101}{39}$', right_node)]
        # offset the nodes
        for node in nodes:
            node.y += offset
        return get_tree_text(nodes, edges)

    def get_right_subtree_text(self, offset):
        """
        @param offset: vertical offset
        @return: a multiline string of TikZ code
        """
        # get the right subtree
        right_nodes, right_edges, right_root = self.get_right_stuff()
        # define the pruned subtree
        left_node = Node('A', -2, 0)
        nodes = right_nodes + [left_node]
        # define the new edge
        edges = right_edges + [(left_node, '$\\frac{52}{21}$', right_root)]
        # offset the nodes
        for node in nodes:
            node.y += offset
        return get_tree_text(nodes, edges)

    def get_tikz_contents(self, show_full_tree, show_pruned_trees):
        out = StringIO()
        if show_full_tree:
            print >> out, self.get_full_tree_text()
            if show_pruned_trees:
                print >> out, self.get_right_subtree_text(-6)
                print >> out, self.get_left_subtree_text(-11)
        elif show_pruned_trees:
            print >> out, self.get_left_subtree_text(0)
            print >> out, self.get_right_subtree_text(-3)
        return out.getvalue().strip()


class SevenLeafLayout(Layout):

    def get_left_stuff(self):
        """
        Use a manual layout.
        @return: (left_nodes, left_edges, left_subtree_root_node)
        """
        j = Node(10, -1, 0)
        i = Node(9, -3, 2)
        h = Node(8, -3, -2)
        a = Node(1, -5, -3)
        b = Node(2, -5, -1)
        c = Node(3, -5, 1)
        d = Node(4, -5, 3)
        left_nodes = [a, b, c, d, h, i, j]
        left_edges = [
                (i, 2, j),
                (d, 2, i),
                (i, 1, c),
                (b, 2, h),
                (j, 3, h),
                (h, 8, a)]
        return left_nodes, left_edges, j

    def get_right_stuff(self):
        """
        Use a manual layout.
        @return: (right_nodes, right_edges, right_subtree_root_node)
        """
        l = Node(12, 2, 0)
        k = Node(11, 4, 2)
        e = Node(5, 6, 3)
        f = Node(6, 6, 1)
        g = Node(7, 4, -1)
        right_nodes = [e, f, g, k, l]
        right_edges = [
                (l, 5, k),
                (k, 5, e),
                (f, 3, k),
                (g, 2, l)]
        return right_nodes, right_edges, l

    def get_full_tree_text(self):
        """
        @return: a multiline string of TikZ code
        """
        # get the left and right subtrees
        left_nodes, left_edges, left_root = self.get_left_stuff()
        right_nodes, right_edges, right_root = self.get_right_stuff()
        # combine the subtrees
        nodes = left_nodes + right_nodes
        edges = left_edges + right_edges + [(left_root, 1, right_root)]
        return get_tree_text(nodes, edges)

    def get_left_subtree_text(self, offset):
        """
        @param offset: vertical offset
        @return: a multiline string of TikZ code
        """
        # get the left subtree
        left_nodes, left_edges, left_root = self.get_left_stuff()
        # define the pruned subtree
        right_node = Node('B', 3, 0)
        nodes = left_nodes + [right_node]
        # define the new edge
        edges = left_edges + [(left_root, '$\\frac{181}{71}$', right_node)]
        # offset the nodes
        for node in nodes:
            node.y += offset
        return get_tree_text(nodes, edges)

    def get_right_subtree_text(self, offset):
        """
        @param offset: vertical offset
        @return: a multiline string of TikZ code
        """
        # get the right subtree
        right_nodes, right_edges, right_root = self.get_right_stuff()
        # define the pruned subtree
        left_node = Node('A', -2, 0)
        nodes = right_nodes + [left_node]
        # define the new edge
        edges = right_edges + [(left_node, '$\\frac{293}{109}$', right_root)]
        # offset the nodes
        for node in nodes:
            node.y += offset
        return get_tree_text(nodes, edges)

    def get_tikz_contents(self, show_full_tree, show_pruned_trees):
        out = StringIO()
        if show_full_tree:
            print >> out, self.get_full_tree_text()
            if show_pruned_trees:
                print >> out, self.get_right_subtree_text(-6)
                print >> out, self.get_left_subtree_text(-11)
        elif show_pruned_trees:
            print >> out, self.get_left_subtree_text(0)
            print >> out, self.get_right_subtree_text(-4)
        return out.getvalue().strip()


def get_tree_text(nodes, edges):
    """
    @param nodes: a collection of node objects
    @param edges: a collection of triples defining edges
    @return: a multiline string of TikZ code
    """
    lines = []
    for node in nodes:
        lines.append(node.get_node_definition())
    for first_node, symbol, second_node in edges:
        edge_object = Edge(symbol, first_node, second_node)
        lines.append(edge_object.get_edge_definition())
    return '\n'.join(lines)


class Edge:

    def __init__(self, symbol, first_node, second_node):
        """
        @param symbol: latex code for the edge label
        @param first_node: an endpoint
        @param second_node: an endpoint
        """
        self.symbol = str(symbol)
        self.first_node = first_node
        self.second_node = second_node

    def get_edge_definition(self):
        """
        @return: a single line defining a TikZ node
        """
        vars = (self.first_node.get_handle(), self.symbol, self.second_node.get_handle())
        line = '\\path (%s) edge node {%s} (%s);' % vars
        return line


class Node:

    def __init__(self, symbol, x, y):
        """
        @param symbol: latex code for the contents of the circle
        @param x: the x position of the node
        @param y: the y position of the node
        """
        self.x = x
        self.y = y
        self.symbol = str(symbol)

    def get_handle(self):
        """
        @return: a unique handle defined by the id of the node
        """
        return 'n' + str(id(self))

    def get_node_definition(self):
        """
        @return: a single line defining a TikZ node
        """
        vars = (self.get_handle(), str(self.x), str(self.y), self.symbol)
        line = '\\node (%s)[draw,shape=circle] at (%s, %s) {%s};' % vars
        return line


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.RadioGroup('nleaves_option', 'use this many leaves', [
                RadioItem('six_leaves', '6'),
                RadioItem('seven_leaves', '7', True)]),
            Form.Float('scaling_factor', 'scaling factor',
                1.0, low_exclusive=0),
            Form.RadioGroup('show_options', 'content options', [
                RadioItem('full_tree_only', 'full tree only', False),
                RadioItem('pruned_trees_only', 'pruned trees only', False),
                RadioItem('all_trees', 'all trees', True)]),
            Form.RadioGroup('output_options', 'output format options', [
                RadioItem('show_tikz', 'show TikZ code only', True),
                RadioItem('show_standalone', 'show contents of a LaTeX file'),
                RadioItem('download_standalone', 'download a LaTeX file'),
                RadioItem('show_pdf', 'view the pdf file'),
                RadioItem('download_pdf', 'download the pdf file'),
                RadioItem('show_png', 'view the png file'),
                RadioItem('download_png', 'download the png file')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # decide which hardcoded tree to use
    if fs.six_leaves:
        layout = SixLeafLayout()
    elif fs.seven_leaves:
        layout = SevenLeafLayout()
    # decide which subtrees to show
    show_full_tree = False
    show_pruned_trees = False
    if fs.full_tree_only:
        show_full_tree = True
    if fs.pruned_trees_only:
        show_pruned_trees = True
    if fs.all_trees:
        show_full_tree = True
        show_pruned_trees = True
    # get the texts
    tikz_text = layout.get_tikz_text(fs.scaling_factor, show_full_tree, show_pruned_trees)
    latex_text = layout.get_latex_text(fs.scaling_factor, show_full_tree, show_pruned_trees)
    # decide the output format
    if fs.show_tikz:
        response_text = tikz_text
        response_headers = [('Content-Type', 'text/plain')]
        return response_headers, response_text
    elif fs.show_standalone:
        response_text = latex_text
        response_headers = [('Content-Type', 'text/plain')]
        return response_headers, response_text
    elif fs.download_standalone:
        response_text = latex_text
        filename = 'tikz-tree.tex'
        contentdisposition = 'attachment'
        response_headers = []
        response_headers.append(('Content-Type', 'text/plain'))
        response_headers.append(('Content-Disposition', "%s; filename=%s" % (contentdisposition, filename)))
        return response_headers, response_text
    elif fs.show_pdf:
        response_text = get_pdf_contents(latex_text)
        response_headers = [('Content-Type', 'application/pdf')]
        return response_headers, response_text
    elif fs.download_pdf:
        response_text = get_pdf_contents(latex_text)
        filename = 'tikz-tree.pdf'
        contentdisposition = 'attachment'
        response_headers = []
        response_headers.append(('Content-Type', 'application/pdf'))
        response_headers.append(('Content-Disposition', "%s; filename=%s" % (contentdisposition, filename)))
        return response_headers, response_text
    elif fs.show_png:
        response_text = get_png_contents(latex_text)
        response_headers = [('Content-Type', 'image/png')]
        return response_headers, response_text
    elif fs.download_png:
        response_text = get_png_contents(latex_text)
        filename = 'tikz-tree.png'
        contentdisposition = 'attachment'
        response_headers = []
        response_headers.append(('Content-Type', 'image/png'))
        response_headers.append(('Content-Disposition', "%s; filename=%s" % (contentdisposition, filename)))
        return response_headers, response_text
    else:
        raise HandlingError('invalid output format')

def create_temp_pdf_file(latex_text):
    """
    @param latex_text: contents of a LaTeX file
    @return: the base of the path name of a temporary pdf file (without the .pdf appended!)
    """
    # write a named temporary latex file
    fd, pathname = tempfile.mkstemp(prefix='webtex', dir='/tmp', text=True)
    fout = os.fdopen(fd, 'w+b')
    fout.write(latex_text)
    fout.close()
    # convert the file to a pdf
    args = ['/usr/bin/pdflatex', '-output-directory', '/tmp', '-interaction', 'nonstopmode', pathname]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    return pathname

def get_png_contents(latex_text):
    """
    This involves using temporary files.
    @param latex_text: contents of a LaTeX file
    @return: contents of a png file
    """
    # create the pdf file
    pathname = create_temp_pdf_file(latex_text)
    # create the png file
    input_arg = pathname + '.pdf'
    output_arg = '-sOutputFile=%s.png' % pathname
    args = ['gs', '-dSAFER', '-dBATCH', '-dNOPAUSE', '-sDEVICE=pnggray', output_arg, input_arg]
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    p_output = p.stdout.read()
    p_error = p.stderr.read()
    # read the png file
    png_pathname = pathname + '.png'
    try:
        fin = open(png_pathname, 'rb')
        png_contents = fin.read()
        fin.close()
    except IOError, e:
        raise HandlingError('failed to create a png file')
    return png_contents

def get_pdf_contents(latex_text):
    """
    This involves using temporary files.
    @param latex_text: contents of a LaTeX file
    @return: contents of a pdf file
    """
    # create the pdf file
    pdf_pathname = create_temp_pdf_file(latex_text) + '.pdf'
    # read the pdf file
    try:
        fin = open(pdf_pathname, 'rb')
        pdf_contents = fin.read()
        fin.close()
    except IOError, e:
        raise HandlingError('failed to create a pdf file')
    return pdf_contents

def main(options):
    layout = SevenLeafLayout()
    print layout.get_latex_text(0.5, True, True)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('--nseconds', dest='nseconds', type='int', default=0, help='seconds to run or 0 to run until ctrl-c')
    #parser.add_option('--ntaxa', dest='ntaxa', type='int', default=20, help='number of taxa in each sampled tree topology')
    options, args = parser.parse_args()
    main(options)
