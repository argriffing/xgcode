"""Construct the rooted tree implied by a phylogenetic contrast matrix.

Each column of the contrast matrix corresponds to a contrast,
and each row corresponds to a taxon.
The order of the columns does not matter.
The sign of each column does not matter;
multiplying one or more columns by -1 does not affect the tree.
The reconstructed tree will be rooted at the center
of the branch implied by the contrast matrix.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import FelTree
import NewickIO
import Contrasts
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    tree = NewickIO.parse(
            Contrasts.g_felsenstein_tree_string, FelTree.NewickTree)
    ordered_labels = ('a', 'b', 'c', 'd', 'e')
    C = Contrasts.get_contrast_matrix(tree, ordered_labels)
    # define the form objects
    form_objects = [
            Form.Matrix('contrast_matrix', 'contrast matrix',
                C, Contrasts.assert_contrast_matrix),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels))]
    return form_objects

def get_form_out():
    return FormOut.Newick()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    C = fs.contrast_matrix
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    # validate the input
    if len(C) != len(ordered_labels):
        raise HandlingError('the number of rows in the contrast matrix should match the number of labels')
    # start to prepare the reponse
    out = StringIO()
    reconstructed_tree = Contrasts.contrast_matrix_to_tree(C, ordered_labels)
    newick_string = NewickIO.get_newick_string(reconstructed_tree)
    print >> out, newick_string
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
