""" Attempt to reject an incorrect tree topology.

If the tree topology could not be rejected,
then provide a certificate of non-rejection.
This is a sequence of internal vertex valuation vectors
which satisfies certain interlacing criteria.
"""

from StringIO import StringIO
import unittest

import Newick
import Form
import FormOut


def get_form():
    """
    @return: the body of a form
    """
    # define default tree strings
    true_tree_string = '((a:1, b:2):3, (c:4, d:5):6, e:7);'
    test_tree_string = '((a, b)x, (c, e)y, d)z;'
    # define the form objects
    form_objects = [
            Form.MultiLine('true_tree', 'true tree', true_tree_string),
            Form.MultiLine('test_tree', 'test topology', test_tree_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # Read the newick trees.
    true_tree = Newick.parse(fs.true_tree, Newick.NewickTree)
    test_tree = Newick.parse(fs.test_tree, Newick.NewickTree)
    out = StringIO()
    print >> out, 'hello, this web thingy is not yet finished'
    return out.getvalue()

def is_branch_compat(nsame, ndifferent, ntarget, nbranches):
    """
    @param nsame: number of placed edges without sign change
    @param ndifferent: number of placed edges with sign change
    @param ntarget: target number of edges with sign change
    @param nbranches: the total number of branches in the tree.
    """
    if nsame + ndifferent > nbranches:
        raise ValueError('branch sign change error')
    if ndifferent > ntarget:
        return False
    npotential = nbranches - nsame
    if npotential < ntarget:
        return False
    return True

def rec_internal(
        id_to_adj, id_to_val,
        nsame, ndifferent, ntarget, nbranches,
        internals, depth):
    """
    This is a recursive function.
    Each level corresponds to an internal vertex.
    Each time a +1/-1 is assigned to an internal vertex,
    check that the number of sign changes on edges is correct.
    @param id_to_adj: node id to list of ids of adjacent nodes
    @param id_to_val: node id to valuation
    @param nsame: number of placed edges without sign change
    @param ndifferent: number of placed edges with sign change
    @param ntarget: target number of edges with sign change
    @param nbranches: the total number of branches in the tree.
    @param internals: list of ids of internal nodes
    @param depth: recursion depth starts at zero
    """
    idcur = internals[depth]
    for value in (-1, 1):
        # Check the number of edges where the signs
        # change and where the signs stay the same
        # under the proposed valuation for the current internal vertex.
        nsame_next = nsame
        ndifferent_next = ndifferent
        for adj in id_to_adj[idcur]:
            adj_val = id_to_val[adj]
            if adj_val is not None:
                prod = adj_val * value
                if prod == -1:
                    ndifferent_next += 1
                elif prod == 1:
                    nsame_next += 1
                else:
                    raise ValueError('edge sign error')
        # If the target number of edges with sign changes
        # is compatible with the current known edge change status
        # then we are OK.
        if is_branch_compat(nsame_next, ndifferent_next, ntarget, nbranches):
            id_to_val[idcur] = value
            if depth == len(internals) - 1:
                yield id_to_val
            else:
                for v in rec_internal(
                        id_to_adj, id_to_val,
                        nsame_next, ndifferent_next, ntarget, nbranches,
                        internals, depth+1):
                    yield v
        # Reset the current value to None.
        id_to_val[idcur] = None

def gen_assignments(
        id_to_adj, id_to_val,
        ntarget, nbranches, internals):
    """
    This is the facade for a recursive function.
    """
    # define the parameters for the recursive function
    # create the generator object
    nsame = 0
    ndifferent = 0
    depth = 0
    obj = rec_internal(
            id_to_adj, id_to_val,
            nsame, ndifferent, ntarget, nbranches,
            internals, depth)
    # return the generator object
    return obj

def rec_eigen():
    """
    This is a recursive function.
    Each level corresponds to an eigenvector.
    """
    pass


class TestThis(unittest.TestCase):

    def test_gen_internal_assignments(self):
        id_to_adj = {
                1 : [6],
                2 : [6],
                3 : [8],
                4 : [7],
                5 : [7],
                6 : [1, 2, 8],
                7 : [4, 5, 8],
                8 : [3, 6, 7]}
        id_to_val = {
                1 : 1,
                2 : 1,
                3 : 1,
                4 : 1,
                5 : -1,
                6 : None,
                7 : None,
                8 : None}
        ntarget = 2
        nbranches = 7
        internals = [6, 7, 8]
        # Get the list of observed compatible assignments.
        ds = list(gen_assignments(
                id_to_adj, id_to_val,
                ntarget, nbranches, internals))
        # Define the only true compatible assignment.
        d = {1: 1, 2: 1, 3: 1, 4: 1, 5: -1, 6: 1, 7: -1, 8: 1}
        # Compare the observed and expected assignments.
        self.assertEqual(len(ds), 1)
        self.assertEqual(ds[0], d)



if __name__ == '__main__':
    unittest.main()
