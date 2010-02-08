"""For each tree, reconstruct the topology from a single eigendecomposition.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import FelTree
import NewickIO
import TreeComparison
import MatrixUtil
import iterutils
from Form import CheckItem
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree lines
    default_tree_lines = [
            '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);',
            '(a:1.062, c:0.190, (d:1.080, b:2.30):2.112);',
            '((d:0.083, b:0.614):0.150, e:0.581, (c:1.290, a:0.359):1.070);',
            '((b:1.749, d:0.523):0.107, e:1.703, (a:0.746, c:0.070):4.025);']
    # define the list of form objects
    form_objects = [
            Form.MultiLine('trees', 'newick trees (one tree per line)',
                '\n'.join(default_tree_lines)),
            Form.Float('epsilon', 'non-negligible eigenvalue', '1e-9'),
            Form.CheckGroup('options', 'show these tree sets', [
                CheckItem('show_all',
                    'all reconstructed trees', True),
                CheckItem('show_incomplete',
                    'incompletely resolved reconstructed trees', True),
                CheckItem('show_conflicting',
                    'conflicting reconstructed trees', True),
                CheckItem('show_negligible',
                    'trees with potentially informative but small loadings',
                    True)])]
    return form_objects


class NegligibleError(Exception):
    """
    This error is raised when a negligible loading is encountered.
    """
    pass


class IncompleteError(Exception):
    """
    This error is raised when a tree reconstruction cannot be completed.
    """
    pass


class AnalysisResult:
    """
    Attempt to reconstruct a tree from the eigendecomposition of the doubly centered distance matrix.
    Report what happens when this is done.
    """
    
    def __init__(self, tree, epsilon):
        """
        @param tree: a newick tree in the felsenstein-inspired format
        @param epsilon: determines whether loadings are considered negligible
        """
        # clear some flags that describe events that occur during reconstruction
        self.is_negligible = False
        self.is_incomplete = False
        self.is_conflicting = False
        # define the trees
        self.tree = tree
        self.reconstructed_tree = None
        # set the threshold for loading negligibility
        self.epsilon = epsilon
        # define some arbitrary ordering of tip names
        self.ordered_names = [node.get_name() for node in tree.gen_tips()]
        # get the distance matrix with respect to this ordering
        D = tree.get_distance_matrix(self.ordered_names)
        # get the Gower doubly centered matrix
        G = MatrixUtil.double_centered(numpy.array(D))
        # get the eigendecomposition of the Gower matrix
        eigenvalues, eigenvector_transposes = np.linalg.eigh(G)
        eigenvectors = eigenvector_transposes.T
        self.sorted_eigensystem = list(reversed(list(sorted((abs(w), v) for w, v in zip(eigenvalues, eigenvectors)))))
        # build the tree recursively using the sorted eigensystem
        indices = set(range(len(self.ordered_names)))
        try:
            # try to reconstruct the tree
            root = self._build_tree(indices, 0)
            root.set_branch_length(None)
            output_tree = Newick.NewickTree(root)
            # convert the tree to the FelTree format
            newick_string = NewickIO.get_newick_string(output_tree)
            self.reconstructed_tree = NewickIO.parse(newick_string, FelTree.NewickTree)
        except NegligibleError:
            self.is_negligible = True
        except IncompleteError:
            self.is_incomplete = True
        else:
            # compare the splits defined by the reconstructed tree to splits in the original tree
            expected_partitions = TreeComparison.get_nontrivial_partitions(self.tree)
            observed_partitions = TreeComparison.get_nontrivial_partitions(self.reconstructed_tree)
            invalid_partitions = observed_partitions - expected_partitions
            if invalid_partitions:
                self.is_conflicting = True

    def _build_tree(self, indices, depth):
        """
        @param indices: a set of indices representing taxa in the current subtree
        @param depth: the depth of the current subtree
        @return: the node representing the subtree
        """
        root = Newick.NewickNode()
        if not indices:
            raise ValueError('trying to build a tree from an empty set of indices')
        elif len(indices) == 1:
            index = list(indices)[0]
            root.set_name(self.ordered_names[index])
        else:
            if depth >= len(self.sorted_eigensystem):
                # the ordered eigenvector loading signs were unable to distinguish each taxon
                raise IncompleteError()
            negative_indices = set()
            positive_indices = set()
            negligible_indices = set()
            w, v = self.sorted_eigensystem[depth]
            for i in indices:
                if abs(v[i]) < self.epsilon:
                    negligible_indices.add(i)
                elif v[i] < 0:
                    negative_indices.add(i)
                else:
                    positive_indices.add(i)
            if negligible_indices:
                # eigenvector loadings near zero are degenerate
                raise NegligibleError()
            for next_indices in (negative_indices, positive_indices):
                if next_indices:
                    child = self._build_tree(next_indices, depth+1)
                    child.set_branch_length(1)
                    root.add_child(child)
                    child.set_parent(root)
        return root

    def get_response_lines(self, options):
        """
        Yield lines that form the result of the analysis.
        @param options: a subset of strings specifying what to show
        """
        preamble_lines = []
        error_lines = []
        if 'show_incomplete' in options and self.is_incomplete:
            error_lines.append('the sequential splits defined by the eigenvectors were insufficient to reconstruct the tree')
        if 'show_conflicting' in options and self.is_conflicting:
            error_lines.append('the reconstructed tree has a split that is incompatible with the original tree')
        if 'show_negligible' in options and self.is_negligible:
            error_lines.append('during reconstruction a negligible eigenvector loading was encountered')
        if 'show_all' in options or error_lines:
            preamble_lines.extend(['original tree:', NewickIO.get_newick_string(self.tree)])
            if self.reconstructed_tree:
                preamble_lines.extend(['reconstructed tree:', NewickIO.get_newick_string(self.reconstructed_tree)])
        return preamble_lines + error_lines


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
    # get the threshold for negligibility of an eigenvector loading
    epsilon = fs.epsilon
    if not (0 <= epsilon < 1):
        raise HandlingError('invalid threshold for negligibility')
    # get the set of selected options
    selected_options = fs.options
    # analyze each tree
    results = []
    for tree in trees:
        results.append(AnalysisResult(tree, epsilon))
    # create the response
    out = StringIO()
    for result in results:
        for line in result.get_response_lines(selected_options):
            print >> out, line
        print >> out
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()


