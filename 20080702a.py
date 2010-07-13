"""Sample some unrooted bifurcating trees.

The number of leaves per tree will be truncated to 52
so that they can be represented by a-z,A-Z,0-9.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import NewickIO
import TreeSampler
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ntrees', 'sample this many trees',
                20, low=1, high=100),
            Form.Integer('leafbase', 'base number of leaves per tree',
                3, low=3, high=20),
            Form.Float('leafmean', 'expected number of extra leaves per tree',
                5, low_inclusive=0, high_inclusive=20),
            Form.Float('branchmean', 'expected length of each branch',
                1, low_exclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # get some samples
    for i in range(fs.ntrees):
        tree = TreeSampler.sample_tree(fs.leafbase, fs.leafmean, fs.branchmean)
        # write the tree
        print >> out, NewickIO.get_newick_string(tree)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
