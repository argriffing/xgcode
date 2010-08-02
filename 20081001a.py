"""Visualize a tree with its centered distance matrix eigendecomposition.

Visualize a tree using the eigendecomposition of its
doubly centered distance matrix.
The line x=0 defines the spectral sign split of the tree.
"""

from scipy import linalg
import numpy as np
import cairo

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import CairoUtil
import MatrixUtil
import Clustering
import NewickIO
import FelTree
import const

g_default_string = const.read('20100730q')

def get_form():
    """
    @return: a list of form objects
    """
    tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    formatted_default_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the list of form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', formatted_default_tree_string),
            Form.CheckGroup('processing_options', 'processing options', [
                Form.CheckItem('internal', 'process internal nodes', True)]),
            Form.CheckGroup('visualization_options', 'visualization options', [
                Form.CheckItem('axes', 'draw axes', True),
                Form.CheckItem('connections', 'draw connections', True)]),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('tree', [])

def get_eigenvectors(row_major_matrix):
    """
    This gets a couple of left eigenvectors
    because of the standard format of rate matrices.
    @param row_major_matrix: this is supposed to be a rate matrix
    @return: a pair of eigenvectors
    """
    R = np.array(row_major_matrix)
    w, vl, vr = linalg.eig(R, left=True, right=True)
    eigenvalue_info = list(sorted((abs(x), i) for i, x in enumerate(w)))
    stationary_eigenvector_index = eigenvalue_info[0][1]
    first_axis_eigenvector_index = eigenvalue_info[1][1]
    second_axis_eigenvector_index = eigenvalue_info[2][1]
    va = vl.T[first_axis_eigenvector_index]
    vb = vl.T[second_axis_eigenvector_index]
    return va, vb

def get_rescaled_vector(v):
    """
    @param v: an array or list of floating point values
    @return: a list of floating point values rescaled to be in the range (0, 1)
    """
    width = max(v) - min(v)
    if not width:
        return [.5 for x in v]
    return [(x-min(v)) / width for x in v]

def get_image(row_major_matrix, incidence_matrix, ordered_names,
        width_and_height, image_format, draw_axes, draw_connections):
    """
    @param row_major_matrix: this is supposed to be a rate matrix
    @param incidence_matrix: for drawing connections
    @param ordered_names: the labels corresponding to rows of the matrix
    @param width_and_height: the dimensions of the output image
    @param image_format: like 'svg', 'png', 'ps', 'pdf', et cetera
    @param draw_axes: True if axes should be drawn
    @param draw_connections: True if connections should be drawn
    @return: a string containing the image data
    """
    width, height = width_and_height
    n = len(row_major_matrix)
    # get eigenvectors scaled to [0, 1]
    va, vb = get_eigenvectors(row_major_matrix)
    rescaled_a = get_rescaled_vector(va)
    rescaled_b = get_rescaled_vector(vb)
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # define the border
    border_fraction = .1
    # draw the axes if requested
    if draw_axes:
        # begin drawing
        context.save()
        context.set_source_rgb(.9, .7, .7)
        # draw the y axis
        dx = max(va) - min(va)
        tx = -min(va)/dx
        xzero = (tx * (1 - 2*border_fraction) + border_fraction) * width
        context.move_to(xzero, 0)
        context.line_to(xzero, height)
        context.stroke()
        # draw the x axis
        dy = max(vb) - min(vb)
        ty = -min(vb)/dy
        yzero = (ty * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(0, yzero)
        context.line_to(width, yzero)
        context.stroke()
        # stop drawing
        context.restore()
    # draw the connections if requested
    if draw_connections:
        # begin drawing
        context.save()
        context.set_source_rgb(.8, .8, .8)
        for i in range(n):
            for j in range(n):
                if not (i < j and incidence_matrix[i][j] > 0):
                    break
                x, y = rescaled_a[i], rescaled_b[i]
                nx = (x * (1 - 2*border_fraction) + border_fraction) * width
                ny = (y * (1 - 2*border_fraction) + border_fraction) * height
                context.move_to(nx, ny)
                x, y = rescaled_a[j], rescaled_b[j]
                nx = (x * (1 - 2*border_fraction) + border_fraction) * width
                ny = (y * (1 - 2*border_fraction) + border_fraction) * height
                context.line_to(nx, ny)
                context.stroke()
        # stop drawing
        context.restore()
    # draw a scatter plot of the states using the eigenvectors as axes
    for i, (x, y) in enumerate(zip(rescaled_a, rescaled_b)):
        state_string = ordered_names[i]
        nx = (x * (1 - 2*border_fraction) + border_fraction) * width
        ny = (y * (1 - 2*border_fraction) + border_fraction) * height
        context.move_to(nx, ny)
        context.show_text(state_string)
    # get the image string
    return cairo_helper.get_image_string()

def get_response_content(fs):
    # start writing the response type
    response_headers = []
    # get the processing options
    use_internal_nodes = fs.internal
    # read the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get the ordered ids and ordered names of the nodes in the tree
    ordered_name_id_pairs = []
    for node in tree.preorder():
        # define the name of the node
        name = ''
        if node.is_tip():
            name = node.get_name()
        # possibly add the node
        if use_internal_nodes:
            ordered_name_id_pairs.append((name, id(node)))
        elif node.is_tip():
            ordered_name_id_pairs.append((name, id(node)))
    ordered_ids = [id_ for name, id_ in ordered_name_id_pairs]
    ordered_names = [name for name, id_ in ordered_name_id_pairs]
    #raise HandlingError('debug: ' + str(ordered_names))
    id_to_index = dict((id_, i) for i, id_ in enumerate(ordered_ids))
    # get the incidence matrix for drawing lines
    n = len(ordered_ids)
    incidence_matrix = [[0]*n for i in range(n)]
    if use_internal_nodes:
        for node in tree.preorder():
            for child in node.gen_children():
                parent_id = id_to_index[id(node)]
                child_id = id_to_index[id(child)]
                incidence_matrix[parent_id][child_id] = 1
                incidence_matrix[child_id][parent_id] = 1
    # get the R matrix from the tree; this is -1/2 times the laplacian matrix
    if use_internal_nodes:
        D = tree.get_full_distance_matrix(ordered_ids)
    else:
        D = tree.get_distance_matrix(ordered_names)
    R_matrix = Clustering.get_R_balaji(D)
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        image_size = (640, 480)
        image_string = get_image(R_matrix, incidence_matrix, ordered_names,
                image_size, ext, fs.axes, fs.connections)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
