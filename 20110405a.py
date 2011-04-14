""" Compute the harmonic extensions of the Schur complement eigenvectors.
"""

from StringIO import StringIO
from collections import defaultdict

import Newick
import Form
import FormOut
import Harmonic
import SeekEigenLacing

g_default_tree = '((1:1, 2:0.5)6:1, (3:0.333333333333, 4:0.5)7:1, 5:1)8;'

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', g_default_tree)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_id_to_name(tree):
    """
    Newick hackery.
    @param tree: a newick-like tree where all nodes are named
    @return: a map from id to name
    """
    return dict((id(x), x.name) for x in tree.preorder())

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

def get_response_content(fs):
    # read the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    # get ids
    leaf_ids = [id(x) for x in tree.gen_tips()]
    internal_ids = [id(x) for x in tree.gen_internal_nodes()]
    nleaves = len(leaf_ids)
    ninternal = len(internal_ids)
    # get map from id to name
    id_to_name = get_id_to_name(tree)
    # Reorder the leaf and the internal node ids according to name order.
    leaf_pair = sorted(
            (id_to_name[x], x) for x in leaf_ids)
    internal_pair = sorted(
            (id_to_name[x], x) for x in internal_ids)
    leaf_ids = zip(*leaf_pair)[1]
    internal_ids = zip(*internal_pair)[1]
    # compute the harmonic extensions
    id_to_val_list = []
    eigenvalues = []
    for i in range(nleaves-1):
        index = i+1
        w, d = Harmonic.get_eigenvalue_and_harmonic_valuations(tree, index)
        id_to_val_list.append(d)
        eigenvalues.append(w)
    # write the report
    out = StringIO()
    print >> out, 'valuations:'
    print >> out
    for i, (w, id_to_val) in enumerate(zip(eigenvalues, id_to_val_list)):
        index = i+1
        print >> out, index, '(1 is Fiedler)'
        print >> out, 'eigenvalue:', w
        print >> out, 'leaves:'
        for x in leaf_ids:
            print >> out, id_to_name[x], id_to_val[x]
        print >> out, 'internal vertices:'
        for x in internal_ids:
            print >> out, id_to_name[x], id_to_val[x]
        print >> out
    return out.getvalue()

