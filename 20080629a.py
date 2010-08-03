"""Given a newick tree, determine if a bipartition is in the split set.

Given a newick tree, determine if a bipartition
is a member of the split set of the tree.
"""

from SnippetUtil import HandlingError
import Util
import FelTree
import NewickIO
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree
    tree_string = NewickIO.daylight_example_tree
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    formatted_tree_string = NewickIO.get_narrow_newick_string(tree, 60)
    # define the selected taxa
    selection = list('ABFG')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree',
                formatted_tree_string),
            Form.MultiLine('selection', 'selected taxa',
                '\n'.join(selection))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get a properly formatted newick tree without caring about branch lengths
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # make sure that each leaf has a unique name
    tip_names = [tip.get_name() for tip in tree.gen_tips()]
    if len(tip_names) < 4:
        raise HandlingError('expected the tree to have at least four leaves')
    if any(name is None for name in tip_names):
        raise HandlingError('each leaf must be named')
    if len(set(tip_names)) != len(tip_names):
        raise HandlingError('each leaf name must be unique')
    # get the selected taxa
    selected_taxa = set(Util.get_stripped_lines(fs.selection.splitlines()))
    # assert that the selected names are actually leaf names
    if set(selected_taxa) - set(tip_names):
        msg = 'one or more selected taxa are not leaf names in the tree'
        raise HandlingError(msg)
    # assert that the selection is not degenerate
    n_selection = len(selected_taxa)
    n_complement = len(tip_names) - n_selection
    if n_selection == 0:
        raise HandlingError('degenerate split: no taxa were selected')
    if n_complement == 0:
        raise HandlingError('degenerate split: all taxa were selected')
    if n_selection == 1:
        msg_a = 'degenerate split: a single selected taxon '
        msg_b = 'can always be separated from the others'
        raise HandlingError(msg_a + msg_b)
    if n_complement == 1:
        msg_a = 'degenerate split: a single selected taxon '
        msg_b = 'can always be separated from the others'
        raise HandlingError(msg_a + msg_b)
    # define the response
    tip_selection = [tip for tip in tree.gen_tips()
            if tip.get_name() in selected_taxa]
    if tree.get_split_branch(tip_selection):
        return 'this split is valid and nontrivial'
    else:
        return 'this split is invalid'
