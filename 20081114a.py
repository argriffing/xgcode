"""
Compare spectral sign splits with and without internal nodes.

For each of a set of trees,
compare spectral sign splits with and without internal nodes.
"""

from StringIO import StringIO
import argparse

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import iterutils
import FelTree
import NewickIO
import TreeComparison
import TreeSampler
import combobreaker

#FIXME use const data

g_cmdline_message = """\
Does the Fiedler split of a certain graph correspond to a cut
of a single branch on the tree from which the graph was derived?
To look for a counterexample to this conjecture,
I am sampling a bunch of random trees.
Force the program to stop running by using Ctrl-c when you get bored."""


def get_form():
    """
    @return: a list of form objects
    """
    # define the default list of arbitrary trees
    default_tree_strings = [
            '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);',
            '(a:1.062, c:0.190, (d:1.080, b:2.30):2.112);',
            '((d:0.083, b:0.614):0.150, e:0.581, (c:1.290, a:0.359):1.070);',
            '((b:1.749, d:0.523):0.107, e:1.703, (a:0.746, c:0.070):4.025);',]
    # define the form objects
    form_objects = [
            Form.MultiLine('trees', 'newick trees (one tree per line)',
                '\n'.join(default_tree_strings))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_vector(distance_matrix):
    """
    Get the principal eigenvector of the double centered distance matrix.
    @param distance_matrix: a distance matrix
    @return: the principal eigenvector of the double centered distance matrix
    """
    D = np.array(distance_matrix)
    n = len(D)
    H = np.eye(n,n) - np.ones((n,n))/n
    HDH = np.dot(H, np.dot(D, H))
    eigenvalues, eigenvector_transposes = np.linalg.eigh(HDH)
    eigenvectors = eigenvector_transposes.T
    max_w, max_v = max([(abs(w), v)
        for w, v in zip(eigenvalues, eigenvectors)])
    return max_v

def get_response_content(fs):
    # get the newick trees.
    trees = []
    for tree_string in iterutils.stripped_lines(fs.trees.splitlines()):
        # parse each tree and make sure that it conforms to various requirements
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        tip_names = [tip.get_name() for tip in tree.gen_tips()]
        if len(tip_names) < 4:
            raise HandlingError(
                    'expected at least four tips '
                    'but found ' + str(len(tip_names)))
        if any(name is None for name in tip_names):
            raise HandlingError('each terminal node must be labeled')
        if len(set(tip_names)) != len(tip_names):
            raise HandlingError('each terminal node label must be unique')
        trees.append(tree)
    # create the response
    out = StringIO()
    same_count = 0
    diff_count = 0
    for tree in trees:
        # make the local paragraph that will be shown if there is an event
        local_out = StringIO()
        has_event = False
        # print the tree
        print >> local_out, NewickIO.get_newick_string(tree)
        # get the tip nodes and the internal nodes
        tip_nodes = []
        internal_nodes = []
        for node in tree.preorder():
            if node.is_tip():
                tip_nodes.append(node)
            else:
                internal_nodes.append(node)
        all_nodes = tip_nodes + internal_nodes
        # get all tip name partitions implied by the tree topology
        valid_partitions = TreeComparison.get_partitions(tree)
        # get results from the augmented distance matrix
        D_full = tree.get_partial_distance_matrix(
                [id(node) for node in all_nodes])
        y_full = get_vector(D_full).tolist()
        y = y_full[:len(tip_nodes)]
        name_selection = frozenset(node.get_name()
                for node, elem in zip(tip_nodes, y) if elem > 0)
        name_complement = frozenset(node.get_name()
                for node, elem in zip(tip_nodes, y) if elem <= 0)
        name_partition_a = frozenset((name_selection, name_complement))
        if name_partition_a not in valid_partitions:
            print >> local_out, 'augmented distance matrix split fail:',
            print >> local_out, name_partition_a
            has_event = True
        # get results from the not-augmented distance matrix
        D = tree.get_partial_distance_matrix([id(node) for node in tip_nodes])
        y = get_vector(D).tolist()
        name_selection = frozenset(node.get_name()
                for node, elem in zip(tip_nodes, y) if elem > 0)
        name_complement = frozenset(node.get_name()
                for node, elem in zip(tip_nodes, y) if elem <= 0)
        name_partition_b = frozenset((name_selection, name_complement))
        if name_partition_b not in valid_partitions:
            print >> local_out, 'not-augmented distance matrix split fail:',
            print >> local_out, name_partition_b
            has_event = True
        # compare the name partitions
        if name_partition_a == name_partition_b:
            same_count += 1
        else:
            diff_count += 1
            print >> local_out, 'this tree was split differently '
            print >> local_out, 'by the different methods:'
            print >> local_out, 'augmented distance matrix split:',
            print >> local_out, name_partition_a
            print >> local_out, 'not-augmented distance matrix split:',
            print >> local_out, name_partition_b
            has_event = True
        # print a newline between trees
        if has_event:
            print >> out, local_out.getvalue()
    # write the summary
    print >> out, 'for this many trees the same split was found:',
    print >> out, same_count
    print >> out, 'for this many trees different splits were found:',
    print >> out, diff_count
    # write the response
    return out.getvalue()

class TrackingChecker:
    def __init__(self, counterexample_filename):
        self.counterexample_filename = counterexample_filename
        self.ncounterexamples = 0
        self.fout = None
    def __call__(self, tree):
        # get the partitions implied by the tree
        valid_partitions = TreeComparison.get_partitions(tree)
        # Get the partition implied by the Fiedler split
        # of the graph derived from the tree.
        tip_nodes = list(tree.gen_tips())
        D = tree.get_partial_distance_matrix(
                [id(node) for node in tip_nodes])
        y = get_vector(D).tolist()
        name_selection = frozenset(node.get_name()
                for node, elem in zip(tip_nodes, y) if elem > 0)
        name_complement = frozenset(node.get_name()
                for node, elem in zip(tip_nodes, y) if elem <= 0)
        name_partition = frozenset((name_selection, name_complement))
        if name_partition not in valid_partitions:
            msg = '\n'.join([
                'invalid partition found:',
                'tree:', NewickIO.get_newick_string(tree),
                'invalid partition:', name_partition])
            if not self.fout:
                self.fout = open(self.counterexample_filename, 'wt')
            print >> self.fout, msg
            print msg
            self.ncounterexamples += 1
        # do not stop looking, even if a counterexample is found
        return False
    def __str__(self):
        return 'found %d counterexamples' % self.ncounterexamples

def gen_trees():
    n_base_leaves = 4
    n_expected_extra_leaves = 1
    expected_branch_length = 1
    while True:
        yield TreeSampler.sample_tree(
                n_base_leaves, n_expected_extra_leaves, expected_branch_length)

def main():
    parser = argparse.ArgumentParser(description=g_cmdline_message)
    parser.add_argument('--counterexamples', default='counterexamples.out',
            help='write counterexamples to this file')
    args = parser.parse_args()
    checker = TrackingChecker(args.counterexamples)
    print combobreaker.run_checker(checker, gen_trees())

if __name__ == '__main__':
    main()

