"""Make a TikZ tree that shows harmonically extended valuations.
"""

from StringIO import StringIO
import numpy as np
import math
import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import tikz
import SpatialTree
import Newick
import NewickIO
import Ftree
import FtreeIO
import FastDaylightLayout

g_default_yaw = 0
g_default_pitch = 0.2
g_default_eigenvector_index = 2
g_radius_inflation = 1.2
g_xy_scale = 20
g_z_scale = 20
g_height_style = 'dashed'

#g_default_tree = '((1:1, 2:0.5)6:1, (3:0.333333333333, 4:0.5)7:1, 5:1)8;' 
g_default_tree = Newick.daylight_example_tree

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree = Newick.parse(g_default_tree, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # get the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', formatted_tree_string),
            Form.Integer('eigenvector_index',
                'eigenvector index (1 is Fiedler)',
                g_default_eigenvector_index, low=1),
            Form.Float('yaw', 'yaw', g_default_yaw),
            Form.Float('pitch', 'pitch', g_default_pitch),
            #Form.Float('scaling_factor',
            #'tikz scaling factor', 1.0, low_exclusive=0),
            Form.TikzFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_v_to_xyz(yaw, v_to_location, v_to_val):
    """
    @param yaw: an angle; rotate the tree layout around its center
    @param v_to_location: maps vertices to (x, y) locations
    @param v_to_val: maps vertices to valuations
    """
    vertices = sorted(v_to_location)
    # force the locations to be centered at the origin
    points = np.array([np.array(loc) for loc in v_to_location.values()])
    center = np.mean(points, axis=0)
    v_to_xy = dict((v, np.array(l) - center) for v, l in v_to_location.items())
    # rotate the locations around their center
    sin_yaw = math.sin(yaw)
    cos_yaw = math.cos(yaw)
    v_to_xyz = {}
    for v in vertices:
        xa, ya = v_to_xy[v]
        xb = xa*cos_yaw - ya*sin_yaw
        yb = xa*sin_yaw + ya*cos_yaw
        zb = v_to_val[v]
        v_to_xyz[v] = np.array([xb, yb, zb])
    return v_to_xyz

def yz_to_y(y, z, pitch):
    cos_pitch = math.cos(pitch)
    sin_pitch = math.sin(pitch)
    return y*sin_pitch + z*cos_pitch

def get_radii(v_to_xyz, pitch):
    """
    The x and y locations are centered.
    @return: horizontal and vertical screen radii for the ellipse
    """
    # first get the distance to the furthest xy
    r_max = max(math.hypot(x,y) for x,y,z in v_to_xyz.values())
    r = r_max * g_radius_inflation
    h_radius = r
    v_radius = yz_to_y(r, 0, pitch)
    return (h_radius, v_radius)

def v_to_shadow(v):
    return ('%s_shadow' % v)

def add_intersection_vertices(T, B, v_to_xyz):
    """
    This is an in-place modification.
    """
    eps = 1e-8
    next_vertex = max(v_to_xyz) + 1
    old_edges = set(T)
    for u_edge in old_edges:
        a, b = u_edge
        ax, ay, az = v_to_xyz[a]
        bx, by, bz = v_to_xyz[b]
        if az * bz < -eps:
            r = next_vertex
            next_vertex += 1
            t = az / (az - bz)
            d = B[u_edge]
            da = t*d
            db = (1-t)*d
            T.remove(u_edge)
            del B[u_edge]
            ea = frozenset((r, a))
            eb = frozenset((r, b))
            T.add(ea)
            T.add(eb)
            B[ea] = da
            B[eb] = db
            rx = t*bx + (1-t)*ax
            ry = t*by + (1-t)*ay
            rz = t*bz + (1-t)*az
            v_to_xyz[r] = np.array([rx, ry, rz])

def xyz_to_tikz_lines(T, B, pitch, v_to_xyz,
        leaves, internal, intersection_vertices):
    """
    The x and y locations are centered.
    The z locations are the harmonically extended valuations.
    @param T: tree topology
    @param B: branch lengths
    @param pitch: an angle; worm eye vs bird eye view of the tree
    @param v_to_xyz: maps vertices to location and valuations
    """
    tikz_lines = []
    plain_vertices = leaves + internal
    # draw shadow vertices
    for v in plain_vertices:
        x, y, z = v_to_xyz[v]
        #style = 'draw,shape=circle,fill=blue,minimum size=3pt'
        style = 'draw'
        line = '\\node (%s)[%s] at (%.4f, %.4f) {};' % (
                v_to_shadow(v), style, x, yz_to_y(y, 0, pitch))
        tikz_lines.append(line)
    # draw vertices of intersection
    for v in intersection_vertices:
        x, y, z = v_to_xyz[v]
        #style = 'draw,shape=circle,fill=blue,minimum size=3pt'
        style = 'draw'
        line = '\\node (%s)[%s] at (%.4f, %.4f) {};' % (
                v, style, x, yz_to_y(y, 0, pitch))
        tikz_lines.append(line)
    # draw the non-positively valuated vertices
    for v in plain_vertices:
        x, y, z = v_to_xyz[v]
        if z >= 0:
            continue
        if v in leaves:
            style = 'draw,shape=circle,fill=black,minimum size=3pt'
        else:
            style = 'draw,shape=circle'
        line = '\\node (%s)[%s] at (%.4f, %.4f) {};' % (
                v, style, x, yz_to_y(y, z, pitch))
        tikz_lines.append(line)
    # draw the edge segments that have non-positive valuation
    for u_edge in T:
        a, b = u_edge
        eps = 1e-8
        if max(v_to_xyz[a][-1], v_to_xyz[b][-1]) >= eps:
            continue
        line = '\\path (%s) edge node {} (%s);' % (a, b)
        tikz_lines.append(line)
    # draw the height bars to non-positively valuated vertices
    for v in plain_vertices:
        x, y, z = v_to_xyz[v]
        if z >= 0:
            continue
        line = '\\path[%s] (%s) edge node {} (%s);' % (
                g_height_style, v, v_to_shadow(v))
        tikz_lines.append(line)
    # draw the translucent ellipse
    h_radius, v_radius = get_radii(v_to_xyz, pitch)
    ellipse_parts = [
            '\\draw[draw=none,fill=lightgray,fill opacity=0.8] (0, 0) ellipse',
            '(%.4fem and %.4fem);' % (h_radius, v_radius)]
    tikz_lines.append(' '.join(ellipse_parts))
    # draw the positively valuated vertices
    for v in plain_vertices:
        x, y, z = v_to_xyz[v]
        if z < 0:
            continue
        if v in leaves:
            style = 'draw,shape=circle,fill=black,minimum size=3pt'
        else:
            style = 'draw,shape=circle'
        line = '\\node (%s)[%s] at (%.4f, %.4f) {};' % (
                v, style, x, yz_to_y(y, z, pitch))
        tikz_lines.append(line)
    # draw the edge segments that have positive valuation
    for u_edge in T:
        a, b = u_edge
        eps = 1e-8
        if max(v_to_xyz[a][-1], v_to_xyz[b][-1]) < eps:
            continue
        line = '\\path (%s) edge node {} (%s);' % (a, b)
        tikz_lines.append(line)
    # draw the height bars to positively valuated vertices
    for v in plain_vertices:
        x, y, z = v_to_xyz[v]
        if z < 0:
            continue
        line = '\\path[%s] (%s) edge node {} (%s);' % (
                g_height_style, v, v_to_shadow(v))
        tikz_lines.append(line)
    return tikz_lines


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

def get_tikz_text(tikz_body, scaling_factor):
    if scaling_factor != 1:
        sf = ',scale=%s' % scaling_factor
    else:
        sf = ''
    #tikz_header = r'\begin{tikzpicture}[auto%s]' % sf
    tikz_header = r'\begin{tikzpicture}[x=1em,y=1em,inner sep=0pt%s]' % sf
    tikz_footer = r'\end{tikzpicture}'
    return '\n'.join([tikz_header, tikz_body, tikz_footer])

def get_tikz_lines(newick, eigenvector_index, yaw, pitch):
    """
    @param eigenvector_index: 1 is Fiedler
    """
    tree = Newick.parse(newick, SpatialTree.SpatialTree) 
    # change the node names and get the new tree string
    for node in tree.preorder():
        node.name = 'n' + str(id(node))
    newick = NewickIO.get_newick_string(tree)
    # do the layout
    layout = FastDaylightLayout.StraightBranchLayout() 
    layout.do_layout(tree) 
    tree.fit((g_xy_scale, g_xy_scale))
    name_to_location = dict((
        x.name, tree._layout_to_display(x.location)) for x in tree.preorder())
    T, B, N = FtreeIO.newick_to_TBN(newick)
    # get some vertices
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    vertices = leaves + internal
    # get the locations
    v_to_location = dict((v, name_to_location[N[v]]) for v in vertices)
    # get the valuations
    w, V = Ftree.TB_to_harmonic_extension(T, B, leaves, internal)
    index_to_val = V[:, eigenvector_index-1]
    v_to_val = dict(
            (vertices[i], g_z_scale*val) for i, val in enumerate(index_to_val))
    # get the coordinates
    v_to_xyz = get_v_to_xyz(yaw, v_to_location, v_to_val)
    # add intersection vertices
    add_intersection_vertices(T, B, v_to_xyz)
    intersection_vertices = sorted(v for v in v_to_xyz if v not in vertices)
    # get lines of the tikz file
    return xyz_to_tikz_lines(T, B, pitch, v_to_xyz,
            leaves, internal, intersection_vertices)

def get_latex_text(tikz_text):
    latex_header = '\n'.join([
        '\\documentclass{article}',
        '\\usepackage{tikz}',
        '\\begin{document}'])
    latex_body = tikz_text
    latex_footer = r'\end{document}'
    return '\n'.join([latex_header, latex_body, latex_footer])

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: the response
    """
    # get the texts
    tikz_lines = get_tikz_lines(
        fs.tree, fs.eigenvector_index, fs.yaw, fs.pitch)
    tikz_text = get_tikz_text('\n'.join(tikz_lines), 1.0)
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

def main(args):
    tikz_lines = get_tikz_lines(
        args.tree, args.eigenvector_index,
        args.yaw, args.pitch)
    tikz_text = get_tikz_text('\n'.join(tikz_lines), 1.0)
    latex_text = get_latex_text(tikz_text)
    print latex_text

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tree', default=g_default_tree,
            help='newick tree with branch lengths')
    parser.add_argument('--eigenvector_index',
            default=g_default_eigenvector_index, type=int,
            help='eigenvector index (1 is Fiedler)')
    parser.add_argument('--yaw', default=g_default_yaw, type=float),
    parser.add_argument('--pitch', default=g_default_pitch, type=float),
    main(parser.parse_args())

