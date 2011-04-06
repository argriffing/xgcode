""" Look for a certificate of non-rejection of a topological alternative.

If the tree topology could not be rejected,
then provide a certificate of non-rejection.
This is a sequence of internal vertex valuation signs
which satisfies certain interlacing criteria.
These criteria are
an interlacing criterion,
sign harmonicity,
and k-cut.
Both trees should have named leaves,
and the set of names should be the same.
The first input tree should have branch lengths.
The second input tree should have named internal vertices.
"""

from StringIO import StringIO
from collections import defaultdict

import Newick
import Form
import FormOut
import Harmonic
import SeekEigenLacing


def get_form():
    """
    @return: the body of a form
    """
    # define default tree strings
    true_s = '(b:2.622, d:1.013, (e:1.496, (a:2.749, c:0.338):1.474):0.889);'
    test_s = '(e, ((c, a)x, d)y, b)z;'
    # define the form objects
    form_objects = [
            Form.MultiLine('true_tree', 'true tree', true_s),
            Form.MultiLine('test_tree', 'test topology', test_s),
            Form.CheckGroup('search_options', 'search criteria', [
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

def get_flags(fs):
    """
    @return: a set of search criterion flags
    """
    flags = set([])
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
    # Read the newick trees.
    true_tree = Newick.parse(fs.true_tree, Newick.NewickTree)
    test_tree = Newick.parse(fs.test_tree, Newick.NewickTree)
    flags = get_flags(fs)
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
    id_to_val_list = []
    for id_to_full_val in id_to_full_val_list:
        d = {}
        for x in true_tree_leaf_ids:
            value = id_to_full_val[x]
            if abs(value) < 1e-8:
                raise ValueError('the true tree is too symmetric')
            elif value < 0:
                s = -1
            else:
                s = 1
            d[id_map[x]] = s
        for x in test_tree_internal_ids:
            d[x] = None
        id_to_val_list.append(d)
    # Attempt to find a sign assignment.
    id_to_vals = SeekEigenLacing.rec_eigen(
            test_id_to_adj, id_to_val_list, flags)
    # Reorder the leaf and the internal node ids according to name order.
    leaf_pair = sorted(
            (test_id_to_name[x], x) for x in test_tree_leaf_ids)
    internal_pair = sorted(
            (test_id_to_name[x], x) for x in test_tree_internal_ids)
    reordered_leaf_ids = zip(*leaf_pair)[1]
    reordered_internal_ids = zip(*internal_pair)[1]
    # Check for a failure to find a certificate.
    if not id_to_vals:
        return 'no non-rejection certificate was found'
    # Start writing the response.
    out = StringIO()
    print >> out, 'leaf sign valuations:'
    for x in reordered_leaf_ids:
        print >> out, test_id_to_name[x], get_sign_string(id_to_vals[x])
    print >> out
    print >> out, 'vertex sign compatible internal vertex valuations:'
    for x in reordered_internal_ids:
        print >> out, test_id_to_name[x], get_sign_string(id_to_vals[x])
    return out.getvalue()
