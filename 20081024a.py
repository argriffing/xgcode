"""Evaluate a split of taxa for a given tree using the exact criterion for deep splits.
"""

from StringIO import StringIO

import numpy

from SnippetUtil import HandlingError
import Form
import NewickIO
import FelTree
import Util
import Clustering

def get_form():
    """
    @return: a list of form objects
    """
    # define default values
    default_tree_string = '((((F:1, E:1):1, D:1):1, C:1):1, B:1, A:1);'
    default_selection = ('A', 'B', 'C')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths', default_tree_string),
            Form.MultiLine('selection', 'selected taxa', '\n'.join(default_selection))]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get the selected names
    selection = list(Util.stripped_lines(StringIO(fs.selection)))
    selected_name_set = set(selection)
    possible_name_set = set(node.get_name() for node in tree.gen_tips())
    extra_names = selected_name_set - possible_name_set
    if extra_names:
        raise HandlingError('the following selected names are not valid tips: %s' % str(tuple(extra_names)))
    complement_name_set = possible_name_set - selected_name_set
    # assert that neither the set of selected names nor its complement is empty
    if not selected_name_set or not complement_name_set:
        raise HandlingError('the selection is degenerate')
    # define an ordering on the tips
    ordered_names = [node.get_name() for node in tree.gen_tips()]
    # convert the selected names to a Y vector
    Y_as_list = []
    for name in ordered_names:
        if name in selected_name_set:
            value = 1
        else:
            value = -1
        Y_as_list.append(value)
    Y = numpy.array(Y_as_list)
    # get the distance matrix
    D = tree.get_distance_matrix(ordered_names)
    # get the R matrix
    R = Clustering.get_R_balaji(D)
    value = numpy.dot(numpy.dot(Y, R), Y.T)
    # report the results
    out = StringIO()
    print >> out, value
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
