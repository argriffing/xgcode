"""Given a distance matrix, get a set of splits using one of several methods.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Util
import BuildTreeTopology
import MatrixUtil
import Form
import FormOut
import const

g_data = const.read('20100730o')

def get_form():
    """
    @return: the body of a form
    """
    # define the default distance matrix and the ordered labels
    lines = Util.get_stripped_lines(g_data.splitlines())
    D = np.array([[float(x) for x in line] for line in lines])
    ordered_labels = list('xyabcmnp')
    # define the form objects
    form_objects = [
            Form.Matrix('matrix', 'distance matrix',
                D, MatrixUtil.assert_predistance),
            Form.MultiLine('labels', 'ordered labels',
                '\n'.join(ordered_labels)),
            Form.RadioGroup('options', 'distance matrix splitting options', [
                Form.RadioItem('option_a',
                    'nj split with nj update', True),
                Form.RadioItem('option_b',
                    'nj split with laplace update'),
                Form.RadioItem('option_c',
                    'spectral with nj fallback and laplace update'),
                Form.RadioItem('option_d',
                    'spectral with partial fallback and laplace update')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def set_to_string(my_set):
    """
    @param my_set: a sequence of arbitrary elements
    """
    return '{' + ', '.join(sorted(str(el) for el in my_set)) + '}'

def split_to_string(my_split):
    """
    @param my_split: a sequence of two taxon sequences
    """
    return ', '.join(set_to_string(taxa) for taxa in my_split)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    D = fs.matrix
    # read the ordered labels
    ordered_labels = Util.get_stripped_lines(StringIO(fs.labels))
    # validate the input
    if len(D) != len(ordered_labels):
        raise HandlingError('the number of taxon labels should match the number of rows in the distance matrix')
    # get the split and update methods
    if fs.option_a:
        split_function = BuildTreeTopology.split_nj
        update_function = BuildTreeTopology.update_nj
    elif fs.option_b:
        split_function = BuildTreeTopology.split_nj
        update_function = BuildTreeTopology.update_using_laplacian
    elif fs.option_c:
        split_function = BuildTreeTopology.split_using_eigenvector_with_nj_fallback
        update_function = BuildTreeTopology.update_using_laplacian
    elif fs.option_d:
        split_function = BuildTreeTopology.split_using_eigenvector
        update_function = BuildTreeTopology.update_using_laplacian
    # get the splits
    index_splits = BuildTreeTopology.get_splits(D, split_function, update_function)
    # start to prepare the reponse
    out = StringIO()
    for index_split in index_splits:
        taxon_split = [[ordered_labels[i] for i in group] for group in index_split]
        print >> out, split_to_string(taxon_split)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
