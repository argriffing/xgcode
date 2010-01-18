import unittest

import day

class TestDay(unittest.TestCase):

    def test_pre_creation(self):
        """
        We cannot set or get items until the tree has been created.
        """
        tree = day.Day()
        self.assertRaises(TypeError, tree.get_x, 3)
        self.assertRaises(TypeError, tree.set_x)
        self.assertRaises(RuntimeError, tree.get_x)
        self.assertRaises(RuntimeError, tree.set_x, 3)

    def test_simple_creation(self):
        """
        Test creation of a simple tree.
        """
        tree = day.Day()
        self.assertRaises(TypeError, tree.end_node, 3)
        self.assertRaises(TypeError, tree.begin_node)
        self.assertRaises(RuntimeError, tree.end_node)
        tree.begin_node(3)
        tree.end_node()
        self.assertRaises(RuntimeError, tree.end_node)
        self.assertRaises(TypeError, tree.get_node_count, 1)
        self.assertEquals(tree.get_node_count(), 1)

    def test_complex_creation(self):
        """
        Test creation of a more complex tree.
        """
        tree = day.Day()
        tree.begin_node(1)
        tree.set_x(42)
        tree.begin_node(10)
        tree.end_node()
        tree.begin_node(20)
        tree.begin_node(21)
        tree.end_node()
        tree.begin_node(22)
        tree.end_node()
        tree.begin_node(23)
        tree.end_node()
        tree.end_node()
        self.assertEquals(tree.get_x(), 42)
        tree.end_node()

    def test_selection(self):
        """
        Test selection of a node by its id.
        """
        # first create a tree
        tree = day.Day()
        tree.begin_node(1)
        tree.set_x(1001)
        tree.begin_node(11)
        tree.set_x(1011)
        tree.end_node()
        tree.begin_node(12)
        tree.set_x(1012)
        tree.end_node()
        tree.begin_node(13)
        tree.end_node()
        tree.end_node()
        # test some stuff without moving the cursor
        self.assertRaises(TypeError, tree.select_node)
        self.assertRaises(TypeError, tree.select_node, 1, 1)
        self.assertRaises(ValueError, tree.select_node, 43)
        # select the root and verify the x value
        tree.select_node(1)
        self.assertEquals(tree.get_x(), 1001)
        # change the selection and verify the x value
        tree.select_node(12)
        self.assertEquals(tree.get_x(), 1012)

    def test_rooting(self):
        # first create a tree
        # 1
        #   11
        #   12
        #     121
        #     122
        #   13
        tree = day.Day()
        tree.begin_node(1)
        tree.set_x(1001)
        tree.begin_node(11)
        tree.set_x(1011)
        tree.end_node()
        tree.begin_node(12)
        tree.set_x(1012)
        tree.begin_node(121)
        tree.end_node()
        tree.begin_node(122)
        tree.end_node()
        tree.end_node()
        tree.begin_node(13)
        tree.end_node()
        tree.end_node()
        # do a sanity check on the created tree
        self.assertEquals(tree.get_node_count(), 6)
        self.assertEquals(tree.get_root_id(), 1)
        # no node should be selected
        self.assertEquals(tree.get_subtree_count(), 0)
        # assert subtree counts of various nodes
        tree.select_node(1)
        self.assertEquals(tree.get_subtree_count(), 3)
        tree.select_node(12)
        self.assertEquals(tree.get_subtree_count(), 2)
        tree.select_node(11)
        self.assertEquals(tree.get_subtree_count(), 0)
        # reroot at a node and make some assertions
        tree.select_node(12)
        tree.reroot()
        self.assertEquals(tree.get_node_count(), 6)
        self.assertEquals(tree.get_subtree_count(), 3)
        self.assertEquals(tree.get_root_id(), 12)



if __name__ == '__main__':
        suite = unittest.TestLoader().loadTestsFromTestCase(TestDay)
        unittest.TextTestRunner(verbosity=2).run(suite)

