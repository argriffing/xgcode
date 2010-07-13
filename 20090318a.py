"""Look for an invalid spectral split when pairwise leaf distances are squared.

Look for a principal coordinate invalid split
when pairwise leaf distances are squared.
Distances between tips of a tree are in some sense
already squared euclidean distances.
The signs of the principal eigenvector of this double centered distance matrix
define a split that is compatible with the tree.
This has not been shown to be the case when the entries of the distance matrix
are transformed by an arbitrary monotonic function,
so this snippet tries to find such a counterexample.
A particular monotonic function of interest is the squaring function,
because the distances may be erroneously squared
if Gower's principal coordinate analysis is applied naively.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import FelTree
import NewickIO
import TreeComparison
import MatrixUtil
import TreeSampler
import iterutils
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree lines
    default_tree_lines = [
            '(((A:1.218, B:0.376):0.62, C:2.586):0.708, D:0.897, E:1.414);',
            '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);',
            '(a:1.062, c:0.190, (d:1.080, b:2.30):2.112);',
            '((d:0.083, b:0.614):0.150, e:0.581, (c:1.290, a:0.359):1.070);',
            '((b:1.749, d:0.523):0.107, e:1.703, (a:0.746, c:0.070):4.025);']
    # define the list of form objects
    form_objects = [
            Form.MultiLine('trees', 'one newick tree per line',
                '\n'.join(default_tree_lines))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def partition_to_string(part):
    """
    @param part: a frozenset of a pair of frozensets of strings
    @return: a string representation of the input
    """
    internal_strings = ['{' + ', '.join(sorted(name_set)) + '}' for name_set in part]
    return '{' + ', '.join(internal_strings) + '}'

def get_principal_eigenvector(M):
    """
    @param M: a 2d numpy array representing a matrix
    @return: the principal eigenvector of M
    """
    eigenvalues, eigenvector_transposes = np.linalg.eigh(M)
    eigenvectors = eigenvector_transposes.T
    eigensystem = [(abs(w), w, v.tolist()) for w, v in zip(eigenvalues, eigenvectors)]
    sorted_eigensystem = list(reversed(sorted(eigensystem)))
    sorted_abs_eigenvalues, sorted_eigenvalues, sorted_eigenvectors = zip(*sorted_eigensystem)
    principal_eigenvector = sorted_eigenvectors[0]
    return principal_eigenvector

def get_principal_coordinate(D):
    """
    Return the principal eigenvector of the pseudoinverse of the laplacian.
    @param D: a distance matrix
    @return: the principal coordinate
    """
    L_pinv = -0.5 * MatrixUtil.double_centered(D)
    return get_principal_eigenvector(L_pinv)

def get_elementwise_squared_matrix(M):
    """
    @param M: a 2d numpy array representing a matrix
    @return: a tranformation of the input matrix where each element is squared
    """
    Q = M.copy()
    nrows, ncols = Q.shape
    for i in range(nrows):
        for j in range(ncols):
            Q[i][j] *= Q[i][j]
    return Q

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the newick trees.
    trees = []
    for tree_string in iterutils.stripped_lines(StringIO(fs.trees)):
        # parse each tree and make sure that it conforms to various requirements
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        tip_names = [tip.get_name() for tip in tree.gen_tips()]
        if len(tip_names) < 4:
            raise HandlingError('expected at least four tips but found ' + str(len(tip_names)))
        if any(name is None for name in tip_names):
            raise HandlingError('each terminal node must be labeled')
        if len(set(tip_names)) != len(tip_names):
            raise HandlingError('each terminal node label must be unique')
        trees.append(tree)
    # begin the response
    out = StringIO()
    # look at each tree
    nerrors = 0
    ncounterexamples = 0
    for tree in trees:
        # get the set of valid partitions implied by the tree
        valid_parts = TreeComparison.get_partitions(tree)
        ordered_tip_names = [tip.get_name() for tip in tree.gen_tips()]
        # assert that the partition implied by the correct formula is valid
        D = np.array(tree.get_distance_matrix(ordered_tip_names))
        loadings = get_principal_coordinate(D)
        nonneg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v >= 0)
        neg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v < 0)
        part = frozenset([nonneg_leaf_set, neg_leaf_set])
        if part not in valid_parts:
            nerrors += 1
            print >> out, 'error: a partition that was supposed to be valid was found to be invalid'
            print >> out, 'tree:', NewickIO.get_newick_string(tree)
            print >> out, 'invalid partition:', partition_to_string(part)
            print >> out
        # check the validity of the partition implied by the incorrect formula
        Q = get_elementwise_squared_matrix(D)
        loadings = get_principal_coordinate(Q)
        nonneg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v >= 0)
        neg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v < 0)
        part = frozenset([nonneg_leaf_set, neg_leaf_set])
        if part not in valid_parts:
            ncounterexamples += 1
            print >> out, 'found a counterexample!'
            print >> out, 'tree:', NewickIO.get_newick_string(tree)
            print >> out, 'invalid partition:', partition_to_string(part)
            print >> out
    print >> out, 'errors found:', nerrors
    print >> out, 'counterexamples found:', ncounterexamples
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

def main():
    filename = 'counterexamples.out'
    fout = open(filename, 'wt')
    print 'Does monotonically transforming the pairwise leaf distances affect the compatibility'
    print 'of the split found using principal coordinate analysis?'
    print 'I am looking through random trees for a tree that is split incompatibly'
    print 'when distances are squared.'
    print 'Use control-c to stop the program when you get bored.'
    try:
        count = 0
        ncounterexamples = 0
        nerrors = 0
        while True:
            count += 1
            # get a random tree
            n_base_leaves = 4
            n_expected_extra_leaves = 1
            expected_branch_length = 1
            tree = TreeSampler.sample_tree(n_base_leaves, n_expected_extra_leaves, expected_branch_length)
            # get the set of valid partitions implied by the tree
            valid_parts = TreeComparison.get_partitions(tree)
            ordered_tip_names = [tip.get_name() for tip in tree.gen_tips()]
            # assert that the partition implied by the correct formula is valid
            D = np.array(tree.get_distance_matrix(ordered_tip_names))
            loadings = get_principal_coordinate(D)
            nonneg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v >= 0)
            neg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v < 0)
            part = frozenset([nonneg_leaf_set, neg_leaf_set])
            if part not in valid_parts:
                nerrors += 1
                print >> fout, 'error: a partition that was supposed to be valid was found to be invalid'
                print >> fout, 'tree:', NewickIO.get_newick_string(tree)
                print >> fout, 'invalid partition:', partition_to_string(part)
                print >> fout
            # check the validity of the partition implied by the incorrect formula
            Q = get_elementwise_squared_matrix(D)
            loadings = get_principal_coordinate(Q)
            nonneg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v >= 0)
            neg_leaf_set = frozenset(tip for tip, v in zip(ordered_tip_names, loadings) if v < 0)
            part = frozenset([nonneg_leaf_set, neg_leaf_set])
            if part not in valid_parts:
                ncounterexamples += 1
                print >> fout, 'found a counterexample!'
                print >> fout, 'tree:', NewickIO.get_newick_string(tree)
                print >> fout, 'invalid partition:', partition_to_string(part)
                print >> fout
    except KeyboardInterrupt, e:
        print 'trees examined:', count
        print 'errors:', nerrors
        print 'counterexamples:', ncounterexamples
    fout.close()

if __name__ == '__main__':
    main()
