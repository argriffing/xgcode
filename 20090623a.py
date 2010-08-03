"""Given a tree, show the sequence of splits found by neighbor joining.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut
import NewickIO
import FelTree
import BuildTreeTopology
import Xtree
import Dendrogram

def get_form():
    """
    @return: the body of a form
    """
    # Define the default tree string with branch lengths
    # and named internal nodes.
    tree_string = '(a:2, (b:2, c:9)g:4, ((d:1, e:3)i:7, f:2)j:1)h;'
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree with branch lengths',
                formatted_tree_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def label_set_to_string(label_set, label_to_name):
    """
    @param label_set: a set of one or more integers
    @param label_to_name: the names associated with the labels
    """
    if len(label_set) == 1:
        label, = label_set
        return label_to_name[label]
    else:
        names = [label_to_name[label] for label in sorted(label_set)]
        return '{' + ', '.join(names) + '}'

class SplitSaver:
    def __init__(self):
        self.splits = []
    def __call__(self, split):
        self.splits.append(split)

def split_to_line(split, names):
    """
    @param split: a split of name indices
    @param names: ordered names of tips
    """
    sides = []
    for indices in split:
        middle = ', '.join(names[i] for i in sorted(indices))
        s = '{' + middle + '}'
        sides.append(s)
    return '{%s, %s}' % tuple(sides)

def do_it_right(D):
    """
    Do neighbor joining correctly.
    @param D: distance matrix
    @return: a sequence of splits
    """
    # use neighbor joining to build the tree, saving the splits in the order they are made
    split_saver = SplitSaver()
    BuildTreeTopology.get_splits(D, BuildTreeTopology.split_nj, BuildTreeTopology.update_nj, split_saver)
    return split_saver.splits

def do_it_wrong(D):
    """
    Do neighbor joining incorrectly.
    @param D: distance matrix
    @return: a sequence of splits
    """
    labels = range(len(D))
    hierarchy = Dendrogram.get_hierarchy(D, BuildTreeTopology.split_nj, labels)
    return set(Dendrogram.hierarchy_to_nontrivial_splits(hierarchy))

def get_response_content(fs):
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # get the distance matrix
    ordered_tip_names = list(sorted(tip.get_name() for tip in tree.gen_tips()))
    D = np.array(tree.get_distance_matrix(ordered_tip_names))
    # begin the output
    out = StringIO()
    # get the order for the correct method of neighbor joining
    splits = do_it_right(D)
    print >> out, 'doing it right'
    print >> out, '\n'.join(split_to_line(split, ordered_tip_names) for split in splits)
    print >> out
    # get the order for the incorrect method of neighbor joining
    splits = do_it_wrong(D)
    print >> out, 'doing it wrong'
    print >> out, '\n'.join(split_to_line(split, ordered_tip_names) for split in splits)
    # write the response
    return out.getvalue()
