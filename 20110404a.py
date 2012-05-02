"""
Look for tree topology rejection power differences.
"""

from StringIO import StringIO
from collections import defaultdict
import time

import numpy as np

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
            Form.CheckGroup('filter_a', 'filter A', [
                Form.CheckItem('always_reject_a',
                    'always reject', False),
                Form.CheckItem('sign_harmonicity_a',
                    'sign harmonicity', True),
                Form.CheckItem('vertex_interlacing_a',
                    'vertex interlacing', True),
                Form.CheckItem('cut_potential_a',
                    'nodal domain cut potential', True),
                Form.CheckItem('orthant_connectivity_a',
                    'orthant connectivity', False)]),
            Form.CheckGroup('filter_b', 'filter B', [
                Form.CheckItem('always_reject_b',
                    'always reject', False),
                Form.CheckItem('sign_harmonicity_b',
                    'sign harmonicity', True),
                Form.CheckItem('vertex_interlacing_b',
                    'vertex interlacing', True),
                Form.CheckItem('cut_potential_b',
                    'nodal domain cut potential', False),
                Form.CheckItem('orthant_connectivity_b',
                    'orthant connectivity', False)]),
            Form.CheckGroup('search_options', 'search options', [
                Form.CheckItem('halt_on_difference',
                    'stop the search when a difference is found', True)])]
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

def get_id_to_name(tree):
    """
    Newick hackery.
    @param tree: a newick-like tree where all nodes are named
    @return: a map from id to name
    """
    return dict((id(x), x.name) for x in tree.preorder())

def get_true_leaf_id_to_test_leaf_id(true_tree, test_tree):
    # Get maps from leaf id to name.
    true_leaf_id_to_name = dict((id(x), x.name) for x in true_tree.gen_tips())
    test_leaf_id_to_name = dict((id(x), x.name) for x in test_tree.gen_tips())
    if len(true_leaf_id_to_name) != len(set(true_leaf_id_to_name.values())):
        raise ValueError('found nonunique leaf names')
    if len(test_leaf_id_to_name) != len(set(test_leaf_id_to_name.values())):
        raise ValueError('found nonunique leaf names')
    true_leaf_name_set = set(true_leaf_id_to_name.values())
    test_leaf_name_set = set(test_leaf_id_to_name.values())
    if true_leaf_name_set != test_leaf_name_set:
        raise ValueError('leaf name mismatch')
    # Get map from name to test tree id.
    name_to_test_leaf_id = dict((b,a) for a, b in test_leaf_id_to_name.items())
    # Get map from true leaf id to test leaf id.
    true_leaf_id_to_test_leaf_id = dict(
            (true_leaf_id, name_to_test_leaf_id[name]) for true_leaf_id,
            name in true_leaf_id_to_name.items())
    return true_leaf_id_to_test_leaf_id

def get_sign_string(arr):
    return ' '.join('+' if x > 0 else '-' for x in arr)

class CheckTreeError(Exception): pass


class CheckTreeData:

    def __init__(self, flags_a, flags_b):
        # store the search parameters
        self.flags_a = flags_a
        self.flags_b = flags_b
        # initialize the search results
        self.acceptance_matrix = np.zeros((2,2))
        self.nerrors = 0
        self.report = None

    def add_acceptances(self, acceptances):
        m = np.zeros_like(self.acceptance_matrix)
        m[acceptances[0]][acceptances[1]] = 1
        self.acceptance_matrix += m

    def add_error(self, e):
        self.nerrors += 1


def check_tree_pair(true_tree, test_tree, data):
    """
    Do this for every pair of sampled trees.
    The returned array has three zeros and a one.
    The position of the one shows
    whether the test tree was accepted or rejected
    with respect to the true tree
    for the pair of methods.
    @param true_tree: a true Newick tree
    @param test_tree: a test Newick tree
    @param data: keeps some running data about the search
    """
    # Get a list of maps from node id to harmonic extension.
    true_tree_leaf_ids = set(id(x) for x in true_tree.gen_tips())
    nleaves = len(true_tree_leaf_ids)
    id_to_full_val_list = [Harmonic.get_harmonic_valuations(
        true_tree, i) for i in range(1, nleaves)]
    id_map =  get_true_leaf_id_to_test_leaf_id(true_tree, test_tree)
    test_id_to_adj = get_id_to_adj(test_tree)
    test_id_to_name = get_id_to_name(test_tree)
    # Get the list of id to val maps with respect to leaf ids of the test tree.
    test_tree_internal_ids = set(id(x) for x in test_tree.gen_internal_nodes())
    test_tree_leaf_ids = set(id(x) for x in test_tree.gen_tips())
    # Fill out the acceptance pair.
    acceptance_pair = []
    for i, flags in enumerate((data.flags_a, data.flags_b)):
        id_to_val_list = []
        for id_to_full_val in id_to_full_val_list:
            d = {}
            for x in true_tree_leaf_ids:
                value = id_to_full_val[x]
                if abs(value) < 1e-8:
                    raise CheckTreeError('the true tree is too symmetric')
                elif value < 0:
                    s = -1
                else:
                    s = 1
                d[id_map[x]] = s
            for x in test_tree_internal_ids:
                d[x] = None
            id_to_val_list.append(d)
        id_to_vals = SeekEigenLacing.rec_eigen(
                test_id_to_adj, id_to_val_list, flags)
        acceptance_pair.append(1 if id_to_vals else 0)
    data.add_acceptances(acceptance_pair)
    found_example = acceptance_pair[0] != acceptance_pair[1]
    if found_example and (data.report is None):
        s = StringIO()
        print >> s, 'true tree:'
        print >> s, true_tree.get_newick_string()
        print >> s
        print >> s, 'test tree:'
        print >> s, test_tree.get_newick_string()
        data.report = s.getvalue()
    return found_example

def get_flags_a(fs):
    """
    @return: a set of search criterion flags
    """
    flags = set([])
    if fs.always_reject_a:
        flags.add(SeekEigenLacing.ALWAYS_REJECT)
    if fs.sign_harmonicity_a:
        flags.add(SeekEigenLacing.SIGN_HARMONICITY)
    if fs.vertex_interlacing_a:
        flags.add(SeekEigenLacing.VERTEX_INTERLACING)
    if fs.cut_potential_a:
        flags.add(SeekEigenLacing.CUT_POTENTIAL)
    if fs.orthant_connectivity_a:
        flags.add(SeekEigenLacing.ORTHANT_CONNECTIVITY)
    return flags

def get_flags_b(fs):
    """
    @return: a set of search criterion flags
    """
    flags = set([])
    if fs.always_reject_b:
        flags.add(SeekEigenLacing.ALWAYS_REJECT)
    if fs.sign_harmonicity_b:
        flags.add(SeekEigenLacing.SIGN_HARMONICITY)
    if fs.vertex_interlacing_b:
        flags.add(SeekEigenLacing.VERTEX_INTERLACING)
    if fs.cut_potential_b:
        flags.add(SeekEigenLacing.CUT_POTENTIAL)
    if fs.orthant_connectivity_b:
        flags.add(SeekEigenLacing.ORTHANT_CONNECTIVITY)
    return flags

def get_response_content(fs):
    flags_a = get_flags_a(fs)
    flags_b = get_flags_b(fs)
    data = CheckTreeData(flags_a, flags_b)
    nseconds = 5
    tm = time.time()
    while time.time() < tm + nseconds:
        # Sample a pair of Newick trees.
        true_f = TreeSampler.sample_tree(fs.nleaves, 0, 1.0)
        test_f = TreeSampler.sample_tree(fs.nleaves, 0, 1.0)
        true_s = NewickIO.get_newick_string(true_f)
        test_s = NewickIO.get_newick_string(test_f)
        true_tree = Newick.parse(true_s, Newick.NewickTree)
        test_tree = Newick.parse(test_s, Newick.NewickTree)
        # Add the pairwise check to the data borg.
        try:
            found_difference = check_tree_pair(true_tree, test_tree, data)
        except CheckTreeError as e:
            data.add_error(e)
        # Check to see if we should stop early.
        if found_difference and fs.halt_on_difference:
            break
    # make the report
    out = StringIO()
    if data.report:
        print >> out, 'found a difference in rejection power'
        print >> out
        print >> out, data.report
        print >> out
    else:
        print >> out, 'failed to find a difference in rejection power'
        print >> out
    print >> out, 'search summary:'
    m = data.acceptance_matrix
    print >> out, 'A reject, B reject:', m[0, 0]
    print >> out, 'A reject, B accept:', m[0, 1]
    print >> out, 'A accept, B reject:', m[1, 0]
    print >> out, 'A accept, B accept:', m[1, 1]
    print >> out, data.nerrors, 'tree symmetry errors'
    return out.getvalue()
