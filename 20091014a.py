"""Compute Fiedler embeddings with various mass vectors on a tree.
"""


import StringIO
import random
import time

import numpy as np
import argparse

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import NewickIO
import FelTree
import Euclid


g_tree_string = '((1:1, 2:0.5)6:1, 5:1, (3:0.33333333333, 4:0.5)7:1)8;'


def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = [
            Form.MultiLine('newick_tree', 'newick tree', g_tree_string)]
    return form_objects

#def get_weighted_embedding(D, m):

#"""
#@param D: a distance matrix
#@param m: a mass vector
#@return: an embedding
#"""
#nvertices = len(m)
#Y = Euclid.weighted_edm_to_points(D, m, nvertices-1)
#sqrtm = np.sqrt(m)
#M = np.diag(sqrtm)
#U, S, VT = np.linalg.svd(np.dot(M, Y))
#Z = np.dot(Y, VT.T)
#return Z

def get_ordered_ids_and_names(tree):
    """
    @param tree: a newick tree
    """
    tip_name_id_pairs = [(node.get_name(), id(node)) for node in tree.gen_tips()]
    internal_name_id_pairs = [(node.get_name(), id(node)) for node in tree.gen_internal_nodes()]
    tip_names, tip_ids = zip(*list(sorted(tip_name_id_pairs)))
    internal_names, internal_ids = zip(*list(sorted(internal_name_id_pairs)))
    ordered_ids = tip_ids + internal_ids
    ordered_names = tip_names + internal_names
    return ordered_ids, ordered_names

def process(tree_string):
    """
    @param tree_string: a newick string
    @return: a multi-line string that summarizes the results
    """
    np.set_printoptions(linewidth=200)
    out = StringIO.StringIO()
    # build the newick tree from the string
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    # get ordered names and ids
    ordered_ids, ordered_names = get_ordered_ids_and_names(tree)
    # get the distance matrix with ordered indices including all nodes in the tree
    nvertices = len(list(tree.preorder()))
    nleaves = len(list(tree.gen_tips()))
    id_to_index = dict((myid, i) for i, myid in enumerate(ordered_ids))
    D = np.array(tree.get_partial_distance_matrix(ordered_ids))
    # define mass vectors
    m_uniform_unscaled = [1]*nvertices
    m_degenerate_unscaled = [1]*nleaves + [0]*(nvertices-nleaves)
    m_uniform = np.array(m_uniform_unscaled, dtype=float) / sum(m_uniform_unscaled)
    m_degenerate = np.array(m_degenerate_unscaled, dtype=float) / sum(m_degenerate_unscaled)
    # show some of the distance matrices
    print >> out, 'ordered names:'
    print >> out, ordered_names
    print >> out
    print >> out, 'embedded points with mass uniformly distributed among all vertices:'
    print >> out, Euclid.edm_to_weighted_points(D, m_uniform)
    print >> out
    print >> out, 'embedded points with mass uniformly distributed among the leaves:'
    print >> out, Euclid.edm_to_weighted_points(D, m_degenerate)
    print >> out
    # return the response
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the response
    result_string = process(fs.newick_tree)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, result_string

if __name__ == '__main__': 
    print process(g_tree_string)
