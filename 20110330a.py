""" Draw an enumeration of possible second cut combos for an invalid tree.

The Fiedler and second cuts are
fiedler: 1+ 2+ 3- 4- 5-
second : 1+ 2+ 3+ 4+ 5-
where the true underlying tree is
((1:1, 2:0.5)6:1, (3:0.333333333333, 4:0.5)7:1, 5:1)8;
and where the topology of the invalid tree for disproof is
((1,2)6, 3, (4,5)7)8;
"""

from SnippetUtil import HandlingError
import Newick
import SpatialTree
import FastDaylightLayout
import Form
import FormOut
import CairoUtil
import DrawEigenLacing

# Get the fourteen combinations of valuations
# defining the cuts of ((1,2)6, 3, (4,5)7)8;
# The (5) node should always have positive valuation;
# this can be renormalized later.
g_annotation = [
            [-1, 1, -1, 1, 1, 1, 1, 1],
            [-1, 1, 1, -1, 1, 1, 1, 1],
            [-1, 1, 1, 1, -1, 1, 1, 1],
            [-1, 1, 1, -1, -1, 1, -1, 1],
            [-3, 3, -1, -1, -1, 3, -1, -1],

            [1, -1, -1, 1, 1, 1, 1, 1],
            [1, -1, 1, -1, 1, 1, 1, 1],
            [1, -1, 1, 1, -1, 1, 1, 1],
            [1, -1, 1, -1, -1, 1, -1, 1],
            [3, -3, -1, -1, -1, 3, -1, -1],

            [-1, -1, -3, 3, 3, -1, 3, 3],
            [-1, -1, 3, -3, 3, -1, 3, 3],
            [-1, -1, 3, 3, -3, -1, 3, 3],
            [-1, -1, 3, -3, -3, -1, -3, 3]]


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.CheckGroup('check_options', 'output options', [
                Form.CheckItem('draw_background', 'draw background', True),
                Form.CheckItem('draw_vertices', 'draw vertices'),
                Form.CheckItem('draw_labels', 'draw labels', True)]),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('example')

def get_response_content(fs):
    # get a properly formatted newick tree with branch lengths
    tree_string = '((1:1, 2:1)6:1, 3:1, (4:1, 5:1)7:1)8;'
    tree = Newick.parse(tree_string, SpatialTree.SpatialTree)
    # get the fiedler-like vertex valuation
    fiedler = [1, 1, -1, -1, -1, 1, -1, -1]
    # create a node id map
    ids = [None]*8
    for node in tree.preorder():
        index = int(node.name) - 1
        ids[index] = id(node)
    # convert fiedler into a dictionary
    v1 = dict((ids[i], float(fiedler[i])) for i in range(8))
    # convert the annotations into dictionaries
    v2s = [dict((ids[i], float(v[i])) for i in range(8)) for v in g_annotation]
    # do the layout
    try:
        layout = FastDaylightLayout.StraightBranchLayout()
        layout.do_layout(tree)
    except RuntimeError, e:
        pass
    # draw the image
    try:
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        return DrawEigenLacing.get_eg2_image(
                tree, (640, 480), ext,
                v1, v2s,
                fs.draw_background, fs.draw_vertices, fs.draw_labels)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
