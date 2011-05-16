"""Draw some TikZ figures to illustrate an ad hoc tree topology resolution.
"""

from StringIO import StringIO

import Form
import FormOut
import tikz

g_segment_length = 1.0
g_cross_radius = 0.2

g_caption = """\
This figure illustrates an ad-hoc tree topology resolution,
given that the signs of the entries of the
$v^2$ (Fiedler) eigenvector of $L*$
partition the leaves like $\\{1,2\\},\\{3,4,5\\}$
while the signs of the entries of the
$v^3$ eigenvector of $L*$
partition the leaves like $\\{1,2,3,4\\},\\{5\\}$.
Because the harmonically extended Fiedler vector cuts a single edge
separating leaves 1 and 2 from the rest of the leaves, we know
that the tree topology is like (a).
This cuts the tree into \'thick\' and \'thin\' nodal domains,
each of which is cut exactly once by the zeros of the harmonically
extended $v^3$ vector.
Because the $v^3$ signs of 1 and 2 are the same,
the single $v^3$ cut of the \'thick\' domain
must be again on the edge that separates leaves 1 and 2
from the rest of the leaves.
Subfigures (b), (c), (d), and (e)
represent the four possible tree topologies given that leaves
1 and 2 are separated by a single branch from the remaining leaves.
The 15 dotted lines of (b), (c), (d), and (e) show potential
locations for the second zero of $v^3$,
given that the second zero must be somewhere in the \'thin\' domain.
Two of these 15 potential locations are illustrated by subfigures (f) and (g).
Subfigure (f) shows $v^3$ zeros which induce
a $\\{1,2,3\\},\\{4,5\\}$ sign partition of the leaves,
which is incompatible
with the observed $\\{1,2,3,4\\},\\{5\\}$ $v^3$ partition.
On the other hand subfigure (g) induces the observed sign partition.
Furthermore, (g) represents the only such combination
of a tree topology and a second $v^3$ cut which induces this sign partition.
Therefore we can deduce that its topology is the only
topology compatible with the tree
from which the $v^2$ and $v^3$ were jointly calculated.
So given only the leaf partitions induced by the
signs of the first two nonconstant eigenvectors of $L*$
we are in this case able to deduce the topology of the tree.
In general this is not always possible.\
"""

#obsolete text from the caption
"""
Of the remaining possible tree topologies (b), (c), (d), and (e),
only topology (d) allows a $v^3$ cut of the \'thin\' domain
which in conjunction
with the previously deduced $v^3$ thick domain cut
sign-isolates leaf 5 from the other leaves
(each path between leaf 5 and each leaf in \\{1,2,3,4\\}
intersects an odd number of
$v^3$ cuts while each path between two leaves in \\{1,2,3,4\\}
intersects an even number of $v^3$ cuts).
"""

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Tikz()

def get_tikz_text(tikz_body):
    tikz_header = '\\begin{tikzpicture}[scale=0.8]'
    tikz_footer = '\\end{tikzpicture}'
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_latex_text(tikz_text):
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\usepackage{subfig}',
        '\\begin{document}'])
    latex_body = tikz_text
    latex_footer = '\\end{document}'
    return '\n'.join([latex_header, latex_body, latex_footer])

def get_float_figure_lines(subfloat_pictures):
    arr = []
    arr.extend([
        '\\begin{figure}',
        #'\\centering'
        '\\begin{center}'
        ])
    for i, picture_text in enumerate(subfloat_pictures):
        arr.extend([
            '\\subfloat[]{',
            '\\label{fig:subfloat-%s}' % i,
            picture_text,
            '}'])
        # add row breaks
        if i in (0, 4):
            arr.append('\\\\')
    arr.extend([
        '\\end{center}',
        '\\caption {',
        g_caption,
        '}',
        '\\label{fig:subfloats}',
        '\\end{figure}'])
    return arr

def get_tikz_uncrossed_line(pt, direction):
    """
    @param pt: a point
    @param direction: a unit vector
    @return: some lines of tikz text
    """
    ax, ay = pt
    dx, dy = direction
    bx = ax + g_segment_length * dx
    by = ay + g_segment_length * dy
    return ['\\draw[color=gray] (%s,%s) -- (%s,%s);' % (ax, ay, bx, by)]

def get_tikz_crossed_line(pt, direction, solid=False):
    """
    @param pt: a point
    @param direction: a unit vector
    @param solid: True if the cross should be solid
    @return: some lines of tikz text
    """
    ax, ay = pt
    dx, dy = direction
    bx = ax + g_segment_length * dx
    by = ay + g_segment_length * dy
    center_x = ax + (g_segment_length / 2.0) * dx
    center_y = ay + (g_segment_length / 2.0) * dy
    cos_90 = 0
    sin_90 = 1
    cos_n90 = 0
    sin_n90 = -1
    qx = center_x + g_cross_radius * (dx*cos_90 - dy*sin_90)
    qy = center_y + g_cross_radius * (dx*sin_90 + dy*cos_90)
    rx = center_x + g_cross_radius * (dx*cos_n90 - dy*sin_n90)
    ry = center_y + g_cross_radius * (dx*sin_n90 + dy*cos_n90)
    style = '' if solid else '[densely dotted]'
    lines = [
            '\n'.join(get_tikz_uncrossed_line(pt, direction)),
            '\\draw%s (%s,%s) -- (%s,%s);' % (style, qx, qy, rx, ry)]
    return lines

def get_fiedler_tikz_lines():
    """
    This should show the Fiedler cut of a tree.
    """
    lines = [
            '\\node[anchor=south] (1) at (1,1) {1};',
            '\\node[anchor=north] (1) at (1,-1) {2};',
            '\\node[anchor=west] (x) at (2.4,0) {$\\{3,4,5\\}$};',
            '\\draw[color=gray] (1,1) -- (1,0);',
            '\\draw[color=gray] (1,0) -- (1,-1);',
            '\\draw[color=gray] (1,0) -- (2,0);',
            '\\draw[color=gray] (2,0) -- (2.4,0.4);',
            '\\draw[color=gray] (2,0) -- (2.4,-0.4);',
            '\\draw[color=gray] (2.4,0.4) -- (2.4,-0.4);',
            '\\draw (1.5,0.2) -- (1.5,-0.2);']
    return lines

def get_preamble_lines():
    """
    @return: some tikz lines common to several tikz figures
    """
    preamble = [
            '\\node[anchor=south] (1) at (1,1) {1};',
            '\\node[anchor=north] (2) at (1,-1) {2};',
            '\\draw[line width=0.1cm,color=gray] (1,1) -- (1,-1);',
            '\\draw[line width=0.1cm,color=gray] (1,0) -- (1.5,0);',
            '\\draw[color=gray] (1.5,0) -- (2,0);',
            '\\draw (1.2,0.2) -- (1.2,-0.2);']
    return preamble

def get_t1_tikz_lines():
    preamble = get_preamble_lines()
    extra_labels = [
            '\\node[anchor=north] (3) at (2,-1) {3};',
            '\\node[anchor=north] (4) at (3,-1) {4};',
            '\\node[anchor=south] (5) at (3,1) {5};']
    crossed_lines = []
    crossed_lines.extend(get_tikz_crossed_line((2,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((2,0), (1,0)))
    crossed_lines.extend(get_tikz_crossed_line((3,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((3,0), (0,1)))
    return preamble + extra_labels + crossed_lines

def get_t2_tikz_lines():
    preamble = get_preamble_lines()
    extra_labels = [
            '\\node[anchor=north] (4) at (2,-1) {4};',
            '\\node[anchor=north] (3) at (3,-1) {3};',
            '\\node[anchor=south] (5) at (3,1) {5};']
    crossed_lines = []
    crossed_lines.extend(get_tikz_crossed_line((2,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((2,0), (1,0)))
    crossed_lines.extend(get_tikz_crossed_line((3,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((3,0), (0,1)))
    return preamble + extra_labels + crossed_lines

def get_t3_tikz_lines():
    preamble = get_preamble_lines()
    extra_labels = [
            '\\node[anchor=north] (5) at (2,-1) {5};',
            '\\node[anchor=north] (3) at (3,-1) {3};',
            '\\node[anchor=south] (4) at (3,1) {4};']
    crossed_lines = []
    crossed_lines.extend(get_tikz_crossed_line((2,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((2,0), (1,0)))
    crossed_lines.extend(get_tikz_crossed_line((3,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((3,0), (0,1)))
    return preamble + extra_labels + crossed_lines

def get_t4_tikz_lines():
    preamble = get_preamble_lines()
    extra_labels = [
            '\\node[anchor=north] (5) at (2,-1) {5};',
            '\\node[anchor=west] (4) at (3,0) {4};',
            '\\node[anchor=south] (3) at (2,1) {3};']
    crossed_lines = []
    crossed_lines.extend(get_tikz_crossed_line((2,0), (0,-1)))
    crossed_lines.extend(get_tikz_crossed_line((2,0), (1,0)))
    crossed_lines.extend(get_tikz_crossed_line((2,0), (0,1)))
    return preamble + extra_labels + crossed_lines

def get_fail_tikz_lines():
    """
    This is an example of a crossing that fails.
    It shows a single failing cut of t1.
    """
    preamble = get_preamble_lines()
    extra_labels = [
            '\\node[anchor=north] (3) at (2,-1) {3};',
            '\\node[anchor=north] (4) at (3,-1) {4};',
            '\\node[anchor=south] (5) at (3,1) {5};']
    lines = []
    lines.extend(get_tikz_crossed_line((2,0), (0,-1), True))
    lines.extend(get_tikz_uncrossed_line((2,0), (1,0)))
    lines.extend(get_tikz_uncrossed_line((3,0), (0,-1)))
    lines.extend(get_tikz_uncrossed_line((3,0), (0,1)))
    return preamble + extra_labels + lines

def get_success_tikz_lines():
    """
    This is an example of a crossing that succeeds.
    It shows a single successful cut of t3.
    """
    preamble = get_preamble_lines()
    extra_labels = [
            '\\node[anchor=north] (5) at (2,-1) {5};',
            '\\node[anchor=north] (3) at (3,-1) {3};',
            '\\node[anchor=south] (4) at (3,1) {4};']
    lines = []
    lines.extend(get_tikz_uncrossed_line((2,0), (0,-1)))
    lines.extend(get_tikz_crossed_line((2,0), (1,0), True))
    lines.extend(get_tikz_uncrossed_line((3,0), (0,-1)))
    lines.extend(get_tikz_uncrossed_line((3,0), (0,1)))
    return preamble + extra_labels + lines

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_fiedler = get_tikz_text('\n'.join(get_fiedler_tikz_lines()))
    tikz_t1 = get_tikz_text('\n'.join(get_t1_tikz_lines()))
    tikz_t2 = get_tikz_text('\n'.join(get_t2_tikz_lines()))
    tikz_t3 = get_tikz_text('\n'.join(get_t3_tikz_lines()))
    tikz_t4 = get_tikz_text('\n'.join(get_t4_tikz_lines()))
    tikz_fail = get_tikz_text('\n'.join(get_fail_tikz_lines()))
    tikz_success = get_tikz_text('\n'.join(get_success_tikz_lines()))
    subfloat_pictures = [
            tikz_fiedler,
            tikz_t1, tikz_t2, tikz_t3, tikz_t4,
            tikz_fail, tikz_success]
    tikz_text = '\n'.join(subfloat_pictures)
    figure_lines = get_float_figure_lines(subfloat_pictures)
    latex_text = get_latex_text('\n'.join(figure_lines))
    # decide the output format
    if fs.tikz:
        return tikz_text
    elif fs.tex:
        return latex_text
    elif fs.pdf:
        return tikz.get_pdf_contents(latex_text)
    elif fs.png:
        return tikz.get_png_contents(latex_text)

