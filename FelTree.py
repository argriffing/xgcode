"""
Make data structures to efficiently represent a tree.
The tree is represented as branches connecting nodes.
This is analogous to edges connecting vertices in other literature.
These data structures are designed to facilitate various
dynamic programming algorithms.
"""

from optparse import OptionParser
import unittest

import NewickIO


class Tree:
    def __init__(self, root=None):
        """
        @param root: the node that is to be the new root or None
        """
        self.set_root(root)

    def get_root(self):
        return self.root

    def set_root(self, node):
        """
        @param node: the node that is to be the new root or None
        """
        if node and not node.is_root():
            path = list(node.gen_path_to_root())
            for child, parent in zip(path[:-1], path[1:]):
                parent.set_directed_branch_to_parent(
                        parent.get_directed_branch_to(child))
            node.set_directed_branch_to_parent(None)
        self.root = node
    
    def get_height(self):
        """
        @return: the number of steps to the deepest node
        """
        max_depth = 0
        stack = [(self.get_root(), 0)]
        while stack:
            node, depth = stack.pop()
            max_depth = max(depth, max_depth)
            for child in node.gen_children():
                stack.append((child, depth+1))
        return max_depth

    def remove_node(self, node):
        """
        Remove any node but the root, and connect its children to its parent.
        @param node: the node to be removed
        """
        if node.is_root():
            raise ValueError('the root of a tree cannot be removed')
        parent = node.get_parent()
        directed_branches_to_children = list(node.gen_exits(parent))
        # remove the connection from the parent to the selected node
        parent.directed_branches.remove(parent.get_directed_branch_to(node))
        # Connect the children of the selected node
        # to the parent of the selected node.
        for directed_branch_to_child in directed_branches_to_children:
            child = directed_branch_to_child.get_target()
            child.get_directed_branch_to_parent().set_target(parent)
            parent.directed_branches.append(directed_branch_to_child)

    def breadth_first(self):
        """
        Generate nodes by a breadth first traversal.
        """
        root = self.get_root()
        if root:
            return list(root.gen_subtree_breadth_first())

    def preorder(self):
        """
        Generate nodes by a preorder traversal.
        """
        return self.breadth_first()

    def postorder(self):
        """
        Generate nodes by a postorder traversal.
        """
        return reversed(list(self.preorder()))

    def gen_non_root_nodes(self):
        for node in self.preorder():
            if not node.is_root():
                yield node

    def gen_internal_nodes(self):
        for node in self.preorder():
            if not node.is_tip():
                yield node

    def gen_tips(self):
        for node in self.preorder():
            if node.is_tip():
                yield node

    def __str__(self):
        """
        Create a printable multi-line string.
        """
        indent = '  '
        lines = []
        stack = [(self.get_root(), 0)]
        while stack:
            node, depth = stack.pop()
            lines.append((indent*depth) + str(node))
            for child in node.gen_children():
                stack.append((child, depth+1))
        return '\n'.join(lines)

    def get_min_binary_branch_length(self):
        """
        This is for checking the Atteson condition for neighbor joining.
        @return: 0 if the tree is not binary, else the smallest branch length.
        """
        min_blen = 0
        max_nneighbors = 0
        for node in self.preorder():
            dneighbors = list(node.get_directed_branches())
            max_nneighbors = max(max_nneighbors, len(dneighbors))
            for dbranch in dneighbors:
                min_branch_length = min(min_blen, dbranch.get_branch_length())
        # if a vertex has more than three neighbors then return zero
        if max_nneighbors > 3:
            return 0
        # otherwise return the smallest branch length
        return min_blen

    def gen_postorder_exits(self):
        """
        Generate (node, directed_branch) pairs.
        The directed branches that target a leaf will be generated first.
        """
        # initialize the free exits associated with each node
        node_id_to_free_exits = {}
        for node in self.preorder():
            node_id = id(node)
            node_id_to_free_exits[node_id] = list(node.gen_directed_branches())
        # initialize the list of exits leading towards dead ends
        exits = []
        for node in self.preorder():
            free_exits = node_id_to_free_exits[id(node)]
            if len(free_exits) == 1:
                exit_branch = free_exits[0]
                target = exit_branch.get_target()
                enter_branch = target.get_directed_branch_to(node)
                exits.append((target, enter_branch))
        # do the iteration
        while exits:
            next_exits = []
            for node, exit_branch in exits:
                yield (node, exit_branch)
                free_exits = node_id_to_free_exits[id(node)]
                free_exits.remove(exit_branch)
                if len(free_exits) == 1:
                    last_exit_branch = free_exits[0]
                    target = last_exit_branch.get_target()
                    enter_branch = target.get_directed_branch_to(node)
                    next_exits.append((target, enter_branch))
                elif len(free_exits) == 0:
                    for next_branch in node.gen_exits(exit_branch.get_target()):
                        target = next_branch.get_target()
                        enter_branch = target.get_directed_branch_to(node)
                        next_exits.append((target, enter_branch))
            exits = next_exits

    def gen_preorder_exits(self):
        """
        Generate (node, directed_branch) pairs.
        The directed branches that target a leaf will be generated last.
        """
        return reversed(list(self.gen_postorder_exits()))

    def get_split_branch(self, selected_tip_list):
        """
        Return a source node and its directed branch.
        Return a source node and its directed branch leading towards the
        unrooted subtree defined by the selected tips.
        The (node, directed branch) pair separaters the selection
        from its complement.
        @param selected_tip_list: list of selected tips defining a bipartition
        @return: None or a (node, directed_branch) pair
        """
        # Define the set of tip ids
        # that are on the target end of the requested branch.
        id_selection = set(id(tip) for tip in selected_tip_list)
        # Define the set of tip ids
        # that are on the target end of the branch with a given id.
        directed_branch_id_to_set = {}
        # fill the dynamic programming data structure
        for node, exiting_branch in self.gen_postorder_exits():
            id_set = set()
            target = exiting_branch.get_target()
            if target.is_tip():
                id_set.update([id(target)])
            else:
                for next_branch in target.gen_exits(node):
                    id_set.update(directed_branch_id_to_set[id(next_branch)])
            if id_set == id_selection:
                return (node, exiting_branch)
            directed_branch_id_to_set[id(exiting_branch)] = id_set


class UndirectedBranch:
    """
    An undirected branch between nodes on the tree.
    This is where the branch length could be stored.
    The color of a branch could also be stored here.
    """
    def __init__(self):
        pass


class Node:
    """
    A leaf or ancestral node in the tree.
    The taxon name and state could be stored here.
    """
    
    UndirectedBranchFactory = UndirectedBranch

    def __init__(self):
        self.directed_branches = []
        self.directed_branch_to_parent = None

    def gen_directed_branches(self):
        return self.directed_branches

    def add_child(self, child):
        """
        If the node to be added has a dummy parent and associated branch then use this branch.
        @param child: the new node to be added as a child of the current node
        """
        # connect the input node to the current node
        child_directed_branch = child.get_directed_branch_to_parent()
        if child_directed_branch:
            # If the node already has a dummy parent
            # then reconnect the child to the current node.
            child_undirected_branch = child_directed_branch.get_undirected_branch()
            child_parent = child_directed_branch.get_target()
            if child_parent:
                raise ValueError(
                        'adding a child node '
                        'that already has a non-degenerate parent')
            child_directed_branch.set_target(self)
        else:
            # If the child has no dummy parent
            # then connect the child to the current node.
            child_undirected_branch = self.UndirectedBranchFactory()
            child_directed_branch = DirectedBranch(self, child_undirected_branch)
            child.add_directed_branch(child_directed_branch)
            child.set_directed_branch_to_parent(child_directed_branch)
        # connect the current node to the child
        self.add_directed_branch(DirectedBranch(child, child_undirected_branch))

    def add_directed_branch(self, directed_branch):
        self.directed_branches.append(directed_branch)

    def set_undirected_branch(self, branch):
        """
        @param branch: the undirected branch that connects with the parent
        """
        if branch.is_root():
            raise ValueError(
                    'setting the undirected branch of the root '
                    'is an undefined operation')
        directed_branch_to_parent = self.get_directed_branch_to_parent()
        parent = branch_to_parent.get_target()
        directed_branch_from_parent = parent.get_directed_branch_to(self)
        directed_branch_to_parent.set_undirected_branch(branch)
        directed_branch_from_parent.set_undirected_branch(branch)

    def set_directed_branch_to_parent(self, directed_branch):
        """
        @param directed_branch: the directed branch leading to the parent
        """
        self.directed_branch_to_parent = directed_branch

    def has_children(self):
        """
        @return: true if the current node has at least one child
        """
        if self.get_child_count():
            return True
        else:
            return False

    def get_child_count(self):
        """
        @return: the number of child nodes
        """
        return len(list(self.gen_children()))

    def get_neighbor_count(self):
        """
        @return: the number of neighbors
        """
        return len(list(self.gen_directed_branches()))

    def get_directed_branch_to_parent(self):
        return self.directed_branch_to_parent

    def get_parent(self):
        """
        @return: the parent node or None
        """
        directed_branch_to_parent = self.get_directed_branch_to_parent()
        if directed_branch_to_parent:
            return directed_branch_to_parent.get_target()

    def get_directed_branch_to(self, neighbor):
        """
        @param neighbor: a neighbor of the current node
        @return: the directed branch to the neighbor
        """
        for directed_branch in self.gen_directed_branches():
            if directed_branch.get_target() is neighbor:
                return directed_branch

    def gen_exits(self, entrance_node):
        """
        Generate all exiting directed branches not leading to the entrance node.
        @param entrance_node: the node which lead to the current node
        """
        for directed_branch in self.gen_directed_branches():
            if directed_branch.get_target() is not entrance_node:
                yield directed_branch

    def gen_neighbors(self):
        for directed_branch in self.gen_directed_branches():
            yield directed_branch.get_target()

    def gen_children(self):
        for directed_branch in self.gen_directed_branches():
            if directed_branch is not self.get_directed_branch_to_parent():
                yield directed_branch.get_target()

    def gen_path_to_root(self):
        """
        Generate all nodes on the path from the current node to the root.
        The current node and the root are included in this path.
        """
        current = self
        while current:
            yield current
            current = current.get_parent()
    
    def gen_subtree_breadth_first(self):
        """
        Note that this is with respect to the rooted tree.
        Generate nodes by a breadth first traversal.
        """
        nodes = [self]
        while nodes:
            next_nodes = []
            for node in nodes:
                yield node
                next_nodes.extend(node.gen_children())
            nodes = next_nodes

    def gen_subtree_preorder(self):
        """
        Note that this is with respect to the rooted tree.
        """
        return list(self.gen_subtree_breadth_first())

    def is_root(self):
        """
        The root is a node with no parent.
        """
        return (not self.get_directed_branch_to_parent())

    def is_tip(self):
        """
        A tip is a node with one neighbor.
        """
        return (len(list(self.gen_directed_branches())) == 1)

    def __str__(self):
        return 'node %d' % id(self)


class DirectedBranch:
    """
    Directed branches are the only connections between the nodes of the tree.
    This extra layer is useful for storing partial results
    during dynamic programming.
    """
    def __init__(self, node, branch):
        """
        @param node: the target of the directed branch
        @param branch: an undirected branch that has information such as length
        """
        assert Node
        self.target_node = node
        self.undirected_branch = branch

    def set_target(self, node):
        self.target_node = node

    def get_target(self):
        return self.target_node

    def set_undirected_branch(self, branch):
        self.undirected_branch = branch

    def get_undirected_branch(self):
        return self.undirected_branch

    def get_branch_length(self):
        return self.undirected_branch.get_branch_length()


class NewickUndirectedBranch(UndirectedBranch):
    def __init__(self, branch_length=None):
        self.set_branch_length(branch_length)

    def set_branch_length(self, branch_length):
        self.branch_length = branch_length

    def get_branch_length(self):
        return self.branch_length

    def create_halves(self):
        """
        Split the branch in half, returning the resulting pair of branches.
        Each of the pair of undirected branches has half the original length.
        @return: a pair of undirected branches
        """
        new_branch_length = self.get_branch_length() / 2.0
        return (NewickUndirectedBranch(new_branch_length),
                NewickUndirectedBranch(new_branch_length))


class NewickNode(Node):
    """
    A Newick node can be named.
    """

    UndirectedBranchFactory = NewickUndirectedBranch

    def __init__(self):
        Node.__init__(self)
        self.name = None

    def get_branch_length(self):
        """
        The branch associated with a node is the one connecting to its parent.
        @return: a floating point branch length or None
        """
        directed_branch = self.get_directed_branch_to_parent()
        if directed_branch:
            return directed_branch.get_undirected_branch().get_branch_length()
        else:
            return None

    def get_name(self):
        return self.name

    def set_branch_length(self, blen):
        """
        The branch associated with a node is the one connecting to its parent.
        If the current node is the root
        then possibly add a dummy parent and an associated undirected branch.
        @return: a floating point branch length or None
        """
        directed_branch = self.get_directed_branch_to_parent()
        if directed_branch:
            directed_branch.get_undirected_branch().set_branch_length(blen)
        else:
            undirected_branch = self.UndirectedBranchFactory(blen)
            directed_branch = DirectedBranch(None, undirected_branch)
            self.add_directed_branch(directed_branch)
            self.set_directed_branch_to_parent(directed_branch)

    def add_name(self, name):
        if self.name:
            raise ValueError('adding a name to a node that is already named')
        else:
            self.name = name

    def set_parent(self, node):
        """
        @param node: the prospective parent
        """
        if node is not self.get_parent():
            raise ValueError('directly changing the parent node is not allowed')

    def __str__(self):
        if self.get_name() is not None:
            return self.get_name()
        else:
            return 'node %d' % id(self)


class NewickTree(Tree):
    """
    A formalization of a phylogenetic tree.
    A Newick tree understands that nodes can have names
    and that branches can have lengths.
    """

    NodeFactory = NewickNode

    def get_affinity_matrix(self, ordered_ids):
        """
        @param ordered_ids: the requested row order by node id
        @return: a row major affinity matrix
        """
        n = len(ordered_ids)
        A = [[0]*n for i in range(n)]
        id_to_index = dict((my_id, i) for i, my_id in enumerate(ordered_ids))
        nodes = list(self.preorder())
        for node_a in nodes:
            index_a = id_to_index.get(id(node_a), None)
            if index_a is None:
                continue
            for directed_branch in node_a.gen_directed_branches():
                node_b = directed_branch.get_target()
                index_b = id_to_index.get(id(node_b), None)
                if index_b is None:
                    continue
                blen = directed_branch.get_undirected_branch().get_branch_length()
                assert blen
                A[index_a][index_b] = 1.0 / blen
        return A

    def get_covariance_matrix(self, ordered_names):
        """
        Unlike the distance matrix, the covariance matrix depends on the root.
        @param ordered_names: the requested row order by name
        @return: a row major covariance matrix
        """
        # get the ordered list of ids of interest, including the root
        name_to_index = dict((name, i) for i, name in enumerate(ordered_names))
        ordered_ids = [None]*len(ordered_names)
        for node in self.gen_tips():
            index = name_to_index.get(node.get_name(), None)
            if index is not None:
                ordered_ids[index] = id(node)
        ordered_ids.append(id(self.get_root()))
        # get the distance matrix that includes distances to the root
        D = self.get_partial_distance_matrix(ordered_ids)
        # transform the distance matrix in place
        for i in range(len(ordered_names)):
            for j in range(len(ordered_names)):
                i_to_root = D[i][len(ordered_names)]
                j_to_root = D[j][len(ordered_names)]
                D[i][j] = 0.5 * (i_to_root + j_to_root - D[i][j])
        # extract the covariance matrix from the transformed distance matrix
        cov = [row[:-1] for row in D[:-1]]
        return cov

    def get_partial_distance_matrix(self, ordered_ids):
        """
        @return: a row major distance matrix
        @param ordered_ids: the requested row order by node id
        """
        # map the id of each node to its index
        id_to_index = dict((id_, i) for i, id_ in enumerate(ordered_ids))
        # get the number of nodes
        n = len(ordered_ids)
        # for each node get the distance to each other node
        distance_matrix = [[0]*n for i in range(n)]
        for node in self.preorder():
            if id(node) not in id_to_index:
                continue
            row = distance_matrix[id_to_index[id(node)]]
            stack = []
            for directed_branch in node.gen_directed_branches():
                next_target = directed_branch.get_target()
                assert next_target
                stack.append((node, next_target, directed_branch.get_undirected_branch().get_branch_length()))
            while stack:
                source, target, distance = stack.pop()
                if id(target) in id_to_index:
                    row[id_to_index[id(target)]] = distance
                for next_branch in target.gen_exits(source):
                    branch_length = next_branch.get_undirected_branch().get_branch_length()
                    next_target = next_branch.get_target()
                    assert next_target, NewickIO.get_newick_string(self)
                    stack.append((target, next_target, distance + branch_length))
        return distance_matrix

    def get_full_distance_matrix(self, ordered_ids=None):
        """
        @return: a row major distance matrix
        @param ordered_ids: the requested row order by node id
        """
        # map the id of each node to its index
        if ordered_ids:
            id_to_index = dict((id_, i) for i, id_ in enumerate(ordered_ids))
        else:
            id_to_index = dict((id(node), i) for i, node in enumerate(self.preorder()))
        # get the number of nodes
        n = len(list(self.preorder()))
        # for each node get the distance to each other node
        distance_matrix = [[0]*n for i in range(n)]
        for node in self.preorder():
            row = distance_matrix[id_to_index[id(node)]]
            stack = []
            for directed_branch in node.gen_directed_branches():
                next_target = directed_branch.get_target()
                assert next_target
                stack.append((node, next_target, directed_branch.get_undirected_branch().get_branch_length()))
            while stack:
                source, target, distance = stack.pop()
                row[id_to_index[id(target)]] = distance
                for next_branch in target.gen_exits(source):
                    branch_length = next_branch.get_undirected_branch().get_branch_length()
                    next_target = next_branch.get_target()
                    assert next_target, NewickIO.get_newick_string(self)
                    stack.append((target, next_target, distance + branch_length))
        return distance_matrix

    def get_distance_matrix(self, ordered_names=None):
        """
        @param ordered_names: the requested order of the names
        @return: a row major distance matrix
        """
        # map the id of each tip to its index
        if ordered_names:
            tip_name_to_index = dict((name, i) for i, name in enumerate(ordered_names))
            tip_id_to_index = dict((id(tip), tip_name_to_index[tip.name]) for tip in self.gen_tips())
        else:
            tip_id_to_index = dict((id(tip), i) for i, tip in enumerate(self.gen_tips()))
        # get the number of tips
        n = len(list(self.gen_tips()))
        # for each tip get the distance to each other tip
        distance_matrix = [[0]*n for i in range(n)]
        for tip in self.gen_tips():
            row = distance_matrix[tip_id_to_index[id(tip)]]
            stack = []
            for directed_branch in tip.gen_directed_branches():
                next_target = directed_branch.get_target()
                assert next_target
                stack.append((tip, next_target, directed_branch.get_undirected_branch().get_branch_length()))
            while stack:
                source, target, distance = stack.pop()
                if target.is_tip():
                    row[tip_id_to_index[id(target)]] = distance
                else:
                    for next_branch in target.gen_exits(source):
                        branch_length = next_branch.get_undirected_branch().get_branch_length()
                        next_target = next_branch.get_target()
                        assert next_target, NewickIO.get_newick_string(self)
                        stack.append((target, next_target, distance + branch_length))
        return distance_matrix

    def get_length(self):
        """
        @return: the sum of the branch lengths of all branches in the tree
        """
        total_length = 0
        for node in self.preorder():
            if not node.is_root():
                branch_length = node.get_directed_branch_to_parent().get_undirected_branch().get_branch_length()
                total_length += branch_length
        return total_length

    def split_branches(self):
        """
        Split each branch in half.
        """
        # get a list of all of the non-root nodes before splitting the branches
        non_root_nodes = [node for node in self.preorder() if not node.is_root()]
        # each non-root node defines a branch
        for node in non_root_nodes:
            # Store the parent of the current node before it is lost.
            parent = node.get_parent()
            # Create the two new undirected branches.
            top_branch, bottom_branch = node.get_directed_branch_to_parent().get_undirected_branch().create_halves()
            # These two connections will be redirected to the new node.
            node_to_parent = node.get_directed_branch_to_parent()
            parent_to_node = parent.get_directed_branch_to(node)
            # Create the new node.
            neo = self.NodeFactory()
            # Redirect the old connections to the new node.
            node_to_parent.set_target(neo)
            parent_to_node.set_target(neo)
            # Create the connections away from the new node.
            neo_to_parent = DirectedBranch(parent, top_branch)
            neo_to_node = DirectedBranch(node, bottom_branch)
            neo.add_directed_branch(neo_to_parent)
            neo.add_directed_branch(neo_to_node)
            neo.set_directed_branch_to_parent(neo_to_parent)
            # Set the undirected branches associated with the pre-existing directed branches.
            node_to_parent.set_undirected_branch(bottom_branch)
            parent_to_node.set_undirected_branch(top_branch)


class SpatialNode(NewickNode):
    """
    A spatial node can have coordinates.
    """

    def __init__(self):
        NewickNode.__init__(self)
        self.location = None

    def set_location(self, location):
        """
        @param location: an (x, y) coordinate pair that defines a point
        """
        self.location = location

    def get_location(self):
        """
        @return: an (x, y) coordinate pair that defines a point
        """
        return self.location


class SpatialTree(NewickTree):
    """
    A spatial tree understands that nodes can have coordinates.
    """

    NodeFactory = SpatialNode

    def rotate(self, origin, theta):
        """
        Rotate each node around the origin.
        @param origin: the (x, y) pair that defines the origin of rotation
        @param theta: the number of radians to rotate each point in the tree about the origin
        """
        pass


def create_test_tree():
    """
    Create a tree for debugging.
    """
    root = Node()
    first_layer = [Node() for i in range(4)]
    for node in first_layer:
        root.add_child(node)
    second_layer = [Node() for i in range(3)]
    for node in second_layer:
        parent = first_layer[0]
        parent.add_child(node)
    return Tree(root)


class Test(unittest.TestCase):

    def test_tree_creation(self):
        create_test_tree()

    def test_tree_stringification(self):
        tree = create_test_tree()
        tree_string = str(tree)

    def test_exit_generator(self):
        tree = create_test_tree()
        # get an observed count
        observed_exit_count = 0
        for node in tree.preorder():
            observed_exit_count += len(node.gen_directed_branches())
        # get a generated count
        generated_exit_count = 0
        for node, directed_branch in tree.gen_postorder_exits():
            generated_exit_count += 1
        self.assertEqual(observed_exit_count, generated_exit_count)

    def test_tip_count(self):
        tree = create_test_tree()
        tip_count = len(list(tree.gen_tips()))
        self.assertEqual(tip_count, 6)

    def test_load_from_string(self):
        tree_string = '((a:1, b:2):3, c:4);'
        tree = NewickIO.parse(tree_string, NewickTree)
        tip_count = len(list(tree.gen_tips()))
        self.assertEqual(tip_count, 3)

    def test_branch_splitting(self):
        tree_string = '((a:1, b:2):3, c:4);'
        tree = NewickIO.parse(tree_string, NewickTree)
        # assert some properties before splitting the branches
        self.assertEqual(len(list(tree.preorder())), 5)
        self.assertEqual(tree.get_length(), 10)
        # split each branch in half, tracking which nodes were added
        old_nodes = list(tree.preorder())
        tree.split_branches()
        new_nodes = [node for node in tree.preorder() if node not in old_nodes]
        # assert some properties after splitting the branches
        self.assertEqual(len(list(tree.preorder())), 9)
        self.assertEqual(tree.get_length(), 10)
        for node in new_nodes:
            exit_branches = list(node.gen_directed_branches())
            self.assertEqual(len(exit_branches), 2)
            branch_lengths = [branch.get_undirected_branch().get_branch_length() for branch in exit_branches]
            self.assertEqual(branch_lengths[0], branch_lengths[1])

    def test_distance_matrix(self):
        tree_string = '((a:1, b:2):3, c:4, d:5);'
        names = list('abcd')
        n = len(names)
        tree = NewickIO.parse(tree_string, NewickTree)
        expected_distance_matrix = [
                [0, 3, 8, 9],
                [3, 0, 9, 10],
                [8, 9, 0, 9],
                [9, 10, 9, 0]]
        distance_matrix = tree.get_distance_matrix(names)
        for i in range(n):
            for j in range(n):
                self.assertAlmostEqual(expected_distance_matrix[i][j], distance_matrix[i][j])

    def test_covariance_matrix(self):
        tree_string = '((a:1, b:2):3, c:4, d:5);'
        tree = NewickIO.parse(tree_string, NewickTree)
        names = list('abcd')
        expected_covariance_matrix = [
                [4, 3, 0, 0],
                [3, 5, 0, 0],
                [0, 0, 4, 0],
                [0, 0, 0, 5]]
        observed_covariance_matrix = tree.get_covariance_matrix(names)
        n = len(names)
        for i in range(n):
            for j in range(n):
                self.assertAlmostEqual(expected_covariance_matrix[i][j], observed_covariance_matrix[i][j])
    
    def test_get_split_branch(self):
        # set up the tree
        tree_string = '((a:1, b:2):3, c:4, d:5);'
        tree = NewickIO.parse(tree_string, NewickTree)
        # look for the branch that separates tips named 'a' and 'b' from the rest of the tree
        tip_selection = [tip for tip in tree.gen_tips() if tip.get_name() in ('a', 'b')]
        node, directed_branch = tree.get_split_branch(tip_selection)
        self.assertEqual(directed_branch.get_undirected_branch().get_branch_length(), 3)
        # look for the branch that separates tips named 'a' and 'c' from the rest of the tree
        tip_selection = [tip for tip in tree.gen_tips() if tip.get_name() in ('a', 'c')]
        result = tree.get_split_branch(tip_selection)
        self.assertEqual(result, None)
        # look for the branch that separates all tips from the rest of the tree
        tip_selection = list(tree.gen_tips())
        result = tree.get_split_branch(tip_selection)
        self.assertEqual(result, None)
        # look for the branch that separates no tips from the rest of the tree
        tip_selection = []
        result = tree.get_split_branch(tip_selection)
        self.assertEqual(result, None)
        # look for the branch that separates the single tip named 'd' from the rest of the tree
        tip_selection = [tip for tip in tree.gen_tips() if tip.get_name() == 'd']
        node, directed_branch = tree.get_split_branch(tip_selection)
        self.assertEqual(directed_branch.get_undirected_branch().get_branch_length(), 5)

    def test_remove_node(self):
        # assert some various preconditions of the tree
        tree = create_test_tree()
        self.assertEqual(len(list(tree.preorder())), 8)
        self.assertEqual(len(list(tree.gen_tips())), 6)
        self.assertEqual(tree.get_height(), 2)
        # removing the most prolific child of the root should flatten the tree
        tree = create_test_tree()
        prolific_child = max((node.get_child_count(), node) for node in tree.get_root().gen_children())[1]
        tree.remove_node(prolific_child)
        self.assertEqual(len(list(tree.preorder())), 7)
        self.assertEqual(len(list(tree.gen_tips())), 6)
        self.assertEqual(tree.get_height(), 1)
        # removing an arbitrary grandchild of the root should not affect the height
        tree = create_test_tree()
        prolific_child = max((node.get_child_count(), node) for node in tree.get_root().gen_children())[1]
        grandchild = list(prolific_child.gen_children())[0]
        tree.remove_node(grandchild)
        self.assertEqual(len(list(tree.preorder())), 7)
        self.assertEqual(len(list(tree.gen_tips())), 5)
        self.assertEqual(tree.get_height(), 2)


class TestNewickDerived(NewickIO.TestNewick):
    TreeFactory = NewickTree


def main():
    print create_test_tree()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        tests = []
        tests.append(unittest.TestLoader().loadTestsFromTestCase(Test))
        tests.append(unittest.TestLoader().loadTestsFromTestCase(TestNewickDerived))
        suite = unittest.TestSuite(tests)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()
