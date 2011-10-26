"""Given a newick tree, give it a new root.

The new root location is specified by giving the parent node of the branch,
the child node of the branch,
and the fraction of the distance between the parent and the child
where the root will be added.
"""

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(a:10, (b:10, c:10)x:10)y;'
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', tree_string),
            Form.SingleLine('parent', 'parent node', 'x'),
            Form.SingleLine('child', 'child node', 'c'),
            Form.Float('fraction',
                'the new root location on this branch', 0.25,
                low_inclusive=0, high_inclusive=1)]
    return form_objects

def get_form_out():
    return FormOut.Newick()

def get_response_content(fs):
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the fraction
    fraction = fs.fraction
    # convert the node names to node objects
    try:
        parent = tree.get_unique_node(fs.parent)
        child = tree.get_unique_node(fs.child)
    except Newick.NewickSearchError as e:
        raise HandlingError(e)
    # allow the parent and child nodes to be specified in the reverse order
    if (parent is not child.parent) and (child is parent.parent):
        parent, child = child, parent
        fraction = 1 - fraction
    # verify the relationship between the parent and child nodes
    if parent is not child.parent:
        msg = 'the given parent and child nodes are not adjacent'
        raise HandlingError(msg)
    # determine the new root node, creating a new one if necessary
    if fraction == 0:
        target = parent
    elif fraction == 1:
        target = child
    else:
        target = Newick.NewickNode()
        tree.insert_node(target, parent, child, fraction)
    if target is tree.root:
        raise HandlingError('the new root is the same as the old root')
    if not target.parent:
        raise HandlingError('topology error')
    # reroot
    old = tree.root
    tree.reroot(target)
    # if the old root has a single child then remove the old root
    if len(old.children) == 1:
        tree.remove_node(old)
    # return the response
    return tree.get_newick_string() + '\n'
