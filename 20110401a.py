""" Make sure that no true topology is rejected.
"""

from StringIO import StringIO
from collections import defaultdict
import numpy as np
import time

import Newick
import Form
import FormOut
import Harmonic
import TreeSampler
import NewickIO
import Newick
import SeekEigenLacing

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nleaves', 'use this many leaves', 5, low=4, high=10),
            Form.CheckGroup('search_options', 'search criteria', [
                Form.CheckItem('always_reject',
                    'always reject', False),
                Form.CheckItem('sign_harmonicity',
                    'sign harmonicity', True),
                Form.CheckItem('vertex_interlacing',
                    'vertex interlacing', True),
                Form.CheckItem('cut_potential',
                    'nodal domain cut potential', True),
                Form.CheckItem('orthant_connectivity',
                    'orthant connectivity', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_id_to_adj(tree):
    """
    Newick hackery.
    @param tree: a newick-like tree
    @return: a map from id to adjacent id list
    """
    id_to_adj = defaultdict(list)
    for a, b in tree.gen_bidirected_branches():
        id_to_adj[id(a)].append(id(b))
    return id_to_adj

class CheckTreeError(Exception): pass

def get_flags(fs):
    """
    @return: a set of search criterion flags
    """
    flags = set([])
    if fs.always_reject:
        flags.add(SeekEigenLacing.ALWAYS_REJECT)
    if fs.sign_harmonicity:
        flags.add(SeekEigenLacing.SIGN_HARMONICITY)
    if fs.vertex_interlacing:
        flags.add(SeekEigenLacing.VERTEX_INTERLACING)
    if fs.cut_potential:
        flags.add(SeekEigenLacing.CUT_POTENTIAL)
    if fs.orthant_connectivity:
        flags.add(SeekEigenLacing.ORTHANT_CONNECTIVITY)
    return flags

def get_response_content(fs):
    flags = get_flags(fs)
    nseconds = 5
    tm = time.time()
    rejected_s = None
    nerrors = 0
    nchecked = 0
    while time.time() < tm + nseconds and not rejected_s:
        nchecked += 1
        # Sample a Newick tree.
        true_f = TreeSampler.sample_tree(fs.nleaves, 0, 1.0)
        true_s = NewickIO.get_newick_string(true_f)
        true_tree = Newick.parse(true_s, Newick.NewickTree)
        # Get the leaf and internal vertex ids for the true tree.
        internal_ids = set(id(x) for x in true_tree.gen_internal_nodes())
        leaf_ids = set(id(x) for x in true_tree.gen_tips())
        nleaves = len(leaf_ids)
        # Get the harmonic valuations for all vertices of the tree.
        id_to_full_val_list = [Harmonic.get_harmonic_valuations(
            true_tree, i) for i in range(1, nleaves)]
        # Check for small valuations at the leaves.
        try:
            for id_to_full_val in id_to_full_val_list:
                for x in leaf_ids:
                    value = id_to_full_val[x]
                    if abs(value) < 1e-8:
                        raise CheckTreeError('the true tree is too symmetric')
        except CheckTreeError as e:
            nerrors += 1
            continue
        # Assign the leaf values and assign None to internal values.
        id_to_val_list = []
        for id_to_full_val in id_to_full_val_list:
            d = {}
            for x in leaf_ids:
                s = -1 if id_to_full_val[x] < 0 else 1
                d[x] = s
            for x in internal_ids:
                d[x] = None
            id_to_val_list.append(d)
        # Define the topology in a different format.
        id_to_adj = get_id_to_adj(true_tree)
        # Check the tree for self-compatibility under the given conditions.
        id_to_list_val = {}
        id_to_vals = SeekEigenLacing.rec_eigen(
            id_to_adj, id_to_val_list, id_to_list_val, 0, flags)
        if not id_to_vals:
            rejected_s = true_s
    # make the report
    out = StringIO()
    if rejected_s:
        print >> out, 'rejected a true tree:'
        print >> out, rejected_s
    else:
        print >> out, 'no true tree was rejected'
    print >> out
    print >> out, nchecked, 'trees were sampled total'
    print >> out, nerrors, 'trees were too symmetric'
    return out.getvalue()
