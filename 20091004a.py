"""Look at the pseudoinverse of a node-weighted Laplacian-like matrix.

Look at the pseudoinverse of a Laplacian-like matrix
for non-uniform node weights.
This Laplacian-like matrix is the cross-product matrix S from the Abdi paper.
"""

from StringIO import StringIO
import random
import time

import numpy as np
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import FelTree
import Euclid
import TreeSampler


def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.Integer('ntaxa', 'number of taxa',
                5, low=3, high=20)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def sample_branch_lengths(tree):
    """
    Modify the tree by setting branch lengths.
    @param tree: a tree
    """
    for node in tree.preorder():
        if not node.is_root():
            branch_length = float(random.randrange(1, 1000))
            node.set_branch_length(branch_length)

def edm_to_S(D, m):
    """
    @param D: a matrix of squared euclidean distances
    @param m: a vector of masses
    @return: an S matrix that is like the pseudoinverse of Laplacian
    """
    n = len(m)
    if D.shape != (n, n):
        raise ValueError('D should be a square matrix conformant to m')
    if any(x < 0 for x in m):
        raise ValueError('each element in m should be nonnegative')
    if not np.allclose(sum(m), 1):
        raise ValueError('the masses should sum to one')
    I = np.eye(n)
    E = I - np.outer(np.ones(n), m)
    S = (-0.5)*np.dot(E, np.dot(D, E))
    return S

def process(ntaxa):
    """
    @param ntaxa: use this many taxa per tree
    @return: a multi-line string that summarizes the results
    """
    np.set_printoptions(linewidth=200)
    # sample an xtree topology
    xtree = TreeSampler.sample_agglomerated_tree(ntaxa)
    # convert the xtree to a FelTree, although I guess this might not be necessary
    tree_string = xtree.get_newick_string()
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    # get ordered ids and the number of leaves and some auxiliary variables
    ordered_ids = get_ordered_ids(tree)
    nleaves = len(list(tree.gen_tips()))
    id_to_index = dict((myid, i) for i, myid in enumerate(ordered_ids))
    # sample random branch lengths
    sample_branch_lengths(tree)
    # get the weighted tree string
    weighted_tree_string = NewickIO.get_newick_string(tree)
    # get the distance matrix relating all vertices
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    # create a mass vector that sums to one
    m = np.array([random.randrange(1, 10) for i in range(len(D))], dtype=float)
    m /= sum(m)
    # get the S matrix
    S = edm_to_S(D, m)
    # get the pseudoinverse of S
    S_pinv = np.linalg.pinv(S)
    # make the response
    out = StringIO()
    print >> out, 'newick tree:', weighted_tree_string
    print >> out
    print >> out, 'm:'
    print >> out, m
    print >> out
    print >> out, 'D:'
    print >> out, D
    print >> out
    print >> out, 'S:'
    print >> out, S
    print >> out
    print >> out, 'pseudoinverse of S:'
    print >> out, S_pinv
    print >> out
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get arguements from the user
    ntaxa = fs.ntaxa
    # get the response
    result_string = process(ntaxa)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, result_string

def get_ordered_ids(tree):
    """
    Maybe I could use postorder here instead.
    @param tree: a tree
    @return: a list of ids beginning with the leaves
    """
    ordered_ids = []
    ordered_ids.extend(id(node) for node in tree.gen_tips())
    ordered_ids.extend(id(node) for node in tree.gen_internal_nodes())
    return ordered_ids

def main(args):
    print process(args.ntaxa)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description=SnippetUtil.docstring_to_title(__doc__))
    parser.add_argument('--ntaxa', type=int, default=5, help='number of taxa in the tree') 
    args = parser.parse_args() 
    main(args) 
