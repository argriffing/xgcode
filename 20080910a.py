"""Evaluate a loss function that compares a query tree to a reference tree.

Note that the distance functions used here are not commutative.
The weighted split distance gives more importance to deep partitions implied by the reference tree.
"""

from SnippetUtil import HandlingError
import NewickIO
import FelTree
import TreeComparison
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    default_tree_string = '((A:1, B:1):1, (C:1, D:1):1, (E:1, F:1):1);'
    # define the form objects
    form_objects = [
            Form.MultiLine('query', 'query tree', default_tree_string),
            Form.MultiLine('reference', 'reference tree', default_tree_string),
            Form.RadioGroup('loss', 'loss function', [
                Form.RadioItem('uniform', 'split distance'),
                Form.RadioItem('weighted', 'weighted split distance', True)]),
            Form.CheckGroup('options', 'normalization options', [
                Form.CheckItem('normalize',
                    'compute the normalized loss', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the query tree
    query_tree = NewickIO.parse(fs.query, FelTree.NewickTree)
    # read the reference tree
    reference_tree = NewickIO.parse(fs.reference, FelTree.NewickTree)
    # calculate the loss using the requested loss function
    if fs.uniform:
        loss_numerator = TreeComparison.get_split_distance(
                query_tree, reference_tree)
    elif fs.weighted:
        loss_numerator = TreeComparison.get_weighted_split_distance(
                query_tree, reference_tree)
    # do the normalization if requested
    if fs.normalize:
        if fs.uniform:
            loss_denominator = float(
                    TreeComparison.get_nontrivial_split_count(reference_tree))
        elif fs.weighted:
            loss_denominator = float(
                    TreeComparison.get_weighted_split_count(reference_tree))
    else:
        loss_denominator = 1
    # return the response
    if loss_denominator:
        return str(loss_numerator / loss_denominator) + '\n'
    else:
        return 'normalization failed\n'

