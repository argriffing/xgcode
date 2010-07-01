"""
This module draws an ascii tree.
"""

class DrawTreeError(Exception): pass

class DrawTree:
    def __init__(self):
        # define some options
        self.vertical_spacing = 3
        self.horizontal_spacing = 6
        self.force_ultrametric = False
        # If branch lengths are used then
        # when the branch length tuning factor is 1.0
        # the drawing should have reasonable width.
        self.use_branch_lengths = False
        self.branch_length_tuning_factor = 1.0

    def _get_branch_length(self, node):
        """
        Get the length of a branch.
        The returned integer represents the number of characters
        in the ascii representation of the branch
        associated with the node.
        @param node: a node in the tree
        @return: an integer length
        """
        if self.use_branch_lengths:
            return int(node.blen * self.branch_length_scaling_factor)
        else:
            return self.horizontal_spacing

    def _add_node_attributes(self, tree, preorder, postorder):
        if self.use_branch_lengths:
            # temporarily add max tip distances postorder
            for node in postorder:
                if node.children:
                    node.max_tip_distance = max(
                            x.blen + x.max_tip_distance for x in node.children)
                else:
                    node.max_tip_distance = 0
            # the max tip distance is the max tip distance of the root
            max_tip_distance = tree.root.max_tip_distance
            # remove the max tip distances postorder
            for node in postorder:
                del node.max_tip_distance
            # Define a branch length scaling factor
            # that depends on the max tip distance.
            sf = (80.0 / max_tip_distance) * self.branch_length_tuning_factor
            self.branch_length_scaling_factor = sf
        # add heights and widths postorder
        for node in postorder:
            children = node.children
            if children:
                gap_height = (len(children) - 1)*self.vertical_spacing 
                child_height = sum(x.height for x in children)
                node.height = gap_height + child_height
                node.width = 1 + max(
                        self._get_branch_length(x) + x.width for x in children)
            else:
                node.height = 1
                node.width = 1 + len(str(self))
        # add offsets preorder
        tree.root.offset = (0, 0)
        for node in preorder:
            row, col = node.offset
            drow = 0
            for child in node.children:
                child.offset = (
                        row + drow, col + 1 + self._get_branch_length(child))
                drow += child.height + self.vertical_spacing
        # add vtargets postorder
        for node in postorder:
            if node.children:
                node.vtarget = (
                        node.children[0].vtarget + node.children[-1].vtarget)/2
            else:
                row, col = node.offset
                node.vtarget = row

    def draw(self, tree):
        """
        @return: a more complicated ascii visualization of the tree than the indented ascii visualization
        """
        # the tree must have a root
        if not tree.root:
            raise DrawTreeError('the tree has no root')
        # cache preorder and postorder node sequences that will not change during drawing
        preorder = tuple(tree.preorder())
        postorder = tuple(tree.postorder())
        # if branch lengths are used then assert that every non-root node has a branch length
        if self.use_branch_lengths:
            for node in preorder:
                if node is not tree.root:
                    if getattr(node, 'blen', None) is None:
                        raise DrawTreeError('found a non-root node without a valid branch length')
        # add a bunch of attributes to the nodes in the tree
        self._add_node_attributes(tree, preorder, postorder)
        # create a rectangle of text
        width = tree.root.width
        height = tree.root.height
        block = dict(((row, col), ' ') for row in range(height) for col in range(width))
        # if we are forcing ultrametric then find the node with the greatest col offset
        if self.force_ultrametric:
            greatest_col_offset = max(node.offset[1] for node in postorder)
        # modify the rectangle of text by drawing on it
        for node in postorder:
            row, col = node.offset
            if not node.children:
                spacer = ''
                if self.force_ultrametric:
                    spacer = '.' * (greatest_col_offset - col)
                s = spacer + ' ' + str(node)
                for i, c in enumerate(s):
                    block[(row, col+i)] = c
            elif len(node.children) == 1:
                for i in range(1 + self._get_branch_length(node.children[0])):
                    block[(row, col+i)] = '-'
            else:
                # draw a long vertical line
                for target_row in range(node.children[0].vtarget, node.children[-1].vtarget):
                    block[(target_row, col)] = '|'
                # draw local intersections and horizontal lines
                for child in node.children:
                    block[(child.vtarget, col)] = '+'
                    for i in range(self._get_branch_length(child)):
                        block[(child.vtarget, col+1+i)] = '-'
                # draw remote intersections
                for child in node.children:
                    if len(child.children) > 1:
                        block[(child.vtarget, col+1+self._get_branch_length(child))] = '+'
        # convert the block of text to lines of text
        lines = [''.join(block[(row, col)] for col in range(width)) for row in range(height)]
        # return a multi-line string with trailing spaces stripped from each line
        return '\n'.join(line.rstrip() for line in lines)

if __name__ == '__main__':
    import Tree
    # create some nodes linked in a tree topology
    root = Tree.TreeNode()
    first_layer = [Tree.TreeNode() for i in range(4)]
    for node in first_layer:
        root.add_child(node)
    second_layer = [Tree.TreeNode() for i in range(3)]
    for node in second_layer:
        first_layer[0].add_child(node)
    # create a tree from these nodes
    tree = Tree.Tree(root)
    # print the tree
    print DrawTree().draw(tree)

