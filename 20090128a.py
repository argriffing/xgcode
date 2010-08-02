"""Compare cuts of a tree before and after removing a set of vertices.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Util
import NewickIO
import Newick
import DrawTree
import FelTree
import Clustering
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree lines
    tree_lines = [
            '(',
            '((a1:1, a2:1):1, (a3:1, a4:1):1):10,',
            '((b1:1, b2:1):1, (b3:1, b4:1):1):20,',
            '((c1:1, c2:1):1, (c3:1, c4:1):1):30);']
    # define the list of form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree',
                '\n'.join(tree_lines)),
            Form.MultiLine('names', 'remove these leaves',
                '\n'.join(('c2', 'c3', 'c4'))),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('show_newick',
                    'show the newick string for each tree', True),
                Form.CheckItem('show_art',
                    'show the ascii art representation of each tree', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_pruned_tree(tree, names_to_remove):
    """
    @param tree: a Newick tree (not a FelTree)
    @param names_to_remove: a set of names of leaves to remove from the tree
    @return: a FelTree
    """
    # get the list of tip nodes to remove
    nodes_to_remove = [node for node in tree.gen_tips() if node.name in names_to_remove]
    # prune the tree
    for node in nodes_to_remove:
        tree.prune(node)
    # merge segmented branches
    internal_nodes_to_remove = [node for node in tree.preorder() if node.get_child_count() == 1]
    for node in internal_nodes_to_remove:
        tree.remove_node(node)
    # convert the tree to the FelTree format
    newick_string = NewickIO.get_newick_string(tree)
    return NewickIO.parse(newick_string, FelTree.NewickTree)

def get_art(tree):
    """
    @param tree: a FelTree
    @return: a multi-line ascii art
    """
    newick_string = NewickIO.get_newick_string(tree)
    simple_tree = NewickIO.parse(newick_string, Newick.NewickTree)
    drawer = DrawTree.DrawTree() 
    drawer.use_branch_lengths = True 
    drawer.force_ultrametric = False 
    drawer.vertical_spacing = 1 
    drawer.horizontal_spacing = 1 
    return drawer.draw(simple_tree)

def process_tree(tree, tree_name, show_newick, show_art):
    """
    @param tree: a FelTree to be split by each method
    @param tree_name: a description of the tree
    @param show_newick: an output option
    @param show_art: an output option
    @return: a multi-line output string
    """
    out = StringIO()
    # be verbose if requested
    if show_newick:
        print >> out, 'newick representation of %s:' % tree_name
        print >> out, Newick.get_narrow_newick_string(tree, 80) 
    if show_art:
        print >> out, 'ascii art representation of %s:' % tree_name
        print >> out, get_art(tree)
    # cut the tree using each method
    ordered_names = list(sorted(node.get_name() for node in tree.gen_tips()))
    n = len(ordered_names)
    D = tree.get_distance_matrix(ordered_names)
    splitters = (Clustering.StoneExactDMS(), Clustering.StoneSpectralSignDMS())
    splitter_names = ('the +1 / -1 split criterion', 'the fiedler criterion')
    for splitter, splitter_name in zip(splitters, splitter_names):
        small_index_selection = splitter.get_selection(D)
        big_index_selection = set(range(n)) - small_index_selection
        names_a = list(sorted(ordered_names[i] for i in small_index_selection))
        names_b = list(sorted(ordered_names[i] for i in big_index_selection))
        print >> out, 'split inferred by %s:' % splitter_name
        print >> out, '{{%s}, {%s}}' % (', '.join(names_a), ', '.join(names_b))
    # return the string
    return out.getvalue()

def get_response_content(fs):
    # get the set of names
    selection = Util.get_stripped_lines(StringIO(fs.names))
    # get the tree
    tree = NewickIO.parse(fs.tree, FelTree.NewickTree)
    # assert that the name selection is compatible with the tree
    selected_name_set = set(selection)
    possible_name_set = set(node.get_name() for node in tree.gen_tips())
    extra_names = selected_name_set - possible_name_set
    if extra_names:
        msg_a = 'the following selected names '
        msg_b = 'are not valid tips: %s' % str(tuple(extra_names)))
        raise HandlingError(msg_a + msg_b)
    # get the pruned tree
    simple_tree = NewickIO.parse(fs.tree, Newick.NewickTree)
    pruned_tree = get_pruned_tree(simple_tree, selected_name_set)
    # begin writing the result
    out = StringIO()
    trees = (tree, pruned_tree)
    tree_names = ('the original tree', 'the pruned tree')
    for tree, tree_name in zip(trees, tree_names):
        print >> out, 'calculating splits of %s:' % tree_name
        print >> out, process_tree(tree, tree_name, fs.show_newick, fs.show_art)
    # return the response
    return out.getvalue()

