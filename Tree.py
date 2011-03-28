"""
A base class for tree structures.
"""


class Tree:
    def __init__(self, root=None):
        self.root = root
        self.cached_preorder = None
        self.cached_postorder = None
    def get_root(self):
        return self.root
    def set_root(self, root):
        self.root = root
    def cache_traversals(self):
        self.cached_preorder = tuple(self.preorder())
        self.cached_postorder = tuple(self.postorder())
    def invalidate_cached_traversals(self):
        self.cached_preorder = None
        self.cached_postorder = None
    def get_max_depth(self):
        return self.root.get_max_depth(0)
    def get_height(self):
        return self.get_max_depth()
    def gen_directed_branches(self):
        v = [self.get_root()]
        while v:
            nextv = []
            for node in v:
                for child in node.gen_children():
                    yield node, child
                    nextv.append(child)
            v = nextv
    def gen_bidirected_branches(self):
        for a, b, in self.gen_directed_branches():
            yield a, b
            yield b, a
    def breadth_first(self):
        if self.root:
            shell = self.root.children
            while shell:
                next_shell = []
                for node in shell:
                    yield node
                    next_shell.extend(node.children)
                shell = next_shell
    def preorder(self):
        if self.cached_preorder is not None:
            for node in self.cached_preorder:
                yield node
        elif self.root:
            for node in self.root.preorder():
                yield node
    def postorder(self):
        if self.cached_postorder is not None:
            for node in self.cached_postorder:
                yield node
        elif self.root:
            for node in self.root.postorder():
                yield node
    def gen_tips(self):
        for node in self.preorder():
            if not node.children:
                yield node
    def gen_internal_nodes(self):
        for node in self.gen_internal_nodes_preorder():
            yield node
    def gen_internal_nodes_preorder(self):
        for node in self.preorder():
            if node.children:
                yield node
    def gen_internal_nodes_postorder(self):
        for node in self.postorder():
            if node.children:
                yield node
    def gen_non_root_nodes(self):
        for node in self.gen_non_root_nodes_preorder():
            yield node
    def gen_non_root_nodes_preorder(self):
        for node in self.preorder():
            if node is not self.root:
                yield node
    def gen_non_root_nodes_postorder(self):
        for node in self.postorder():
            if node is not self.root:
                yield node
    def get_path_to_root(self, node):
        path = []
        while node:
            path.append(node)
            node = node.parent
        return path
    def add_depths(self):
        if self.root:
            self.root.depth = 0
            for node in self.preorder():
                for child in node.children:
                    child.depth = node.depth + 1
    def del_depths(self):
        for node in self.preorder():
            del node.depth
    def gen_description_lines(self):
        """
        Yield single line strings describing features of the tree.
        """
        self.cache_traversals()
        yield 'the number of nodes is %d' % len(list(self.preorder()))
        yield 'the number of child nodes of the root is %d' % len(self.root.children)
        yield 'the number of nodes with no children is %d' % len(list(self.gen_tips()))
        yield 'the number of nodes with at least one child is %d' % len(list(self.gen_internal_nodes_preorder()))
        yield 'the number of nodes with exactly one child is %d' % len([node for node in self.preorder() if len(node.children) == 1])
        yield 'the number of nodes with more than three children is %d' % len([node for node in self.preorder() if len(node.children) > 3])
        yield 'the number of non-root nodes with exactly three children is %d' % len([node for node in self.preorder() if len(node.children) == 3 and node is not self.root])
        yield 'the maximum depth of the tree is %d' % self.get_max_depth()
    def __str__(self):
        self.add_depths()
        arr = [('  '*node.depth) + str(node) for node in self.preorder()]
        return '\n'.join(arr)

class TreeNode:
    def __init__(self):
        self.parent = None
        self.children = []
    def get_parent(self):
        return self.parent
    def has_children(self):
        return any(self.children)
    def get_neighbor_count(self):
        count = self.get_child_count()
        if self.parent:
            count += 1
        return count
    def get_child_count(self):
        return len(self.get_children())
    def get_children(self):
        return self.children
    def gen_children(self):
        return self.children
    def set_parent(self, parent):
        self.parent = parent
    def add_child(self, child):
        self.children.append(child)
    def gen_tips(self):
        for node in self.preorder():
            if not node.children:
                yield node
    def preorder(self):
        stack = [self]
        while stack:
            node = stack.pop()
            yield node
            stack.extend(node.children)
    def gen_subtree_preorder(self):
        """
        Note that this is with respect to the rooted tree.
        """
        return list(self.preorder())
    def postorder(self):
        for node in reversed(list(self.preorder())):
            yield node
    def get_max_depth(self, depth):
        if not self.children:
            return depth
        return 1 + max(child.get_max_depth(depth) for child in self.children)


if __name__ == '__main__':
    root = TreeNode()
    first_layer = [TreeNode() for i in range(4)]
    for node in first_layer:
        root.add_child(node)
    second_layer = [TreeNode() for i in range(3)]
    for node in second_layer:
        first_layer[0].add_child(node)
    tree = Tree(root)
    print tree

