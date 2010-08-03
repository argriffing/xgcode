"""Seek a counterexample to an eigenvector sign sufficiency conjecture.

Look for a specific counterexample
to the eigenvector sign sufficiency conjecture.
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

#TODO allow the user to define the string and the target sign pattern

g_epsilon = 1e-10

# no known branch lengths for this tree generate the target sign pattern
g_tree_string = '((1, 2), 3, (4, 5));'

# this tree has known branch lengths that generate the target sign pattern
#g_tree_string = '((1, 2), 5, (3, 4));'

g_target_sign_patterns = np.array([
    [-1, -1, +1, +1, +1],
    [-1, -1, -1, -1, +1],
    [+1, -1, +1, +1, +1],
    [-1, -1, +1, -1, -1]])

class CounterexampleError(Exception): pass

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    form_objects = []
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

def get_response_content(fs):
    # on the web we have a short attention span
    return process(nseconds=2) + '\n'

def process(nseconds=None):
    """
    @param nseconds: allow this many seconds to run or None to run forever
    @return: a multi-line string that summarizes the results
    """
    # load the tree
    tree = NewickIO.parse(g_tree_string, FelTree.NewickTree)
    # get the alphabetically ordered tip names
    ordered_tip_names = list(sorted(node.get_name() for node in tree.gen_tips()))
    # initialize the search
    start_time = time.time()
    nsamples_rejected = 0
    nsamples_accepted = 0
    counterexample_message = 'no counterexample was found'
    try:
        while True:
            elapsed_time = time.time() - start_time
            if nseconds and elapsed_time > nseconds:
                break
            # sample some random branch lengths
            sample_branch_lengths(tree)
            # get the distance matrix
            D = np.array(tree.get_distance_matrix(ordered_tip_names))
            # get the projections onto the MDS axes of the leaves
            X = Euclid.edm_to_points(D)
            # if any coordinate is near zero then reject the sample
            if np.min(np.abs(X)) < g_epsilon:
                nsamples_rejected += 1
                continue
            # see if the sign pattern matches for each coordinate
            for v_observed, v_target in zip(X.T, g_target_sign_patterns):
                hadamard_product = v_observed * v_target
                all_positive = all(x>0 for x in hadamard_product)
                all_negative = all(x<0 for x in hadamard_product)
                if not (all_positive or all_negative):
                    # the target sign pattern was not met
                    break
            else:
                # the sign pattern matched for each coordinate so we have a counterexample
                msg = NewickIO.get_newick_string(tree)
                raise CounterexampleError(msg)
            # increment the count of accepted samples
            nsamples_accepted += 1
    except KeyboardInterrupt, e:
        pass
    except CounterexampleError, e:
        counterexample_message = str(e)
    # make the response
    out = StringIO()
    print >> out, elapsed_time, 'seconds elapsed'
    print >> out, 'epsilon for nonzero coordinates:', g_epsilon
    print >> out, 'search status:', counterexample_message
    print >> out, nsamples_rejected, 'samples were rejected because of a near-zero coordinate'
    print >> out, nsamples_accepted, 'samples were accepted'
    return out.getvalue().strip()

def main(args):
    print 'looking for branch lengths that give the target sign pattern...'
    print process(args.nseconds)

if __name__ == '__main__': 
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--nseconds', type=int, default=0, help='seconds to run or 0 to run until ctrl-c') 
    args = parser.parse_args() 
    main(args) 
