"""
Neighborhood joining related algorithms.
This is like neighbor joining,
but it can split more than two branches from the star tree at each iteration.
"""

from StringIO import StringIO
import unittest

import Newick
import Clustering
import NeighborJoining
import MatrixUtil
import Util
import JC69
import NewickIO
import FelTree


class NeighborhoodJoiningError(Exception):
    pass


class SequentialNumberGenerator:
    def __init__(self, start=0):
        self.current = start
    def get_next(self):
        next = self.current
        self.current += 1
        return next


def split_distance_matrix(D, selection, complement):
    """
    @param D: a distance matrix
    @param selection: the selected indices
    @return: a pair of distance matrices
    """
    # get the crossing distances
    v = get_crossing_distances(D, selection, complement)
    # get the output matrices
    output_matrices = []
    for A in (selection, complement):
        rows = []
        for i in sorted(A):
            row = []
            for j in sorted(A):
                row.append(D[i][j])
            row.append(v[i])
            rows.append(row)
        row = []
        for i in sorted(A):
            row.append(v[i])
        row.append(0)
        rows.append(row)
        output_matrices.append(rows)
    return output_matrices

def get_crossing_distances(D, selection, complement):
    """
    This is another method used to get the length of a branch joining two subtrees.
    @param D: a distance matrix
    @param selection: one of the clusters, represented by a set of indices
    @param complement: the other cluser
    @return: the distance from each vertex to the root of the other subtree
    """
    # get the number of vertices
    n = len(selection | complement)
    # get the expected distance to a point in the other cluster
    E = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            E[i] = sum(D[i][j] for j in B)/len(B)
    # get the mean of E
    E_bar = sum(E)/n
    # get the vector of distances to a virtual point that may or may not be on the tree
    v = [0]*n
    for A, B in ((selection, complement), (complement, selection)):
        for i in A:
            v[i] = E[i] - len(A)*E_bar/n
    # get the branch length
    ma = min(v[i] + v[j] - D[i][j] for i in selection for j in selection)/2
    mb = min(v[i] + v[j] - D[i][j] for i in complement for j in complement)/2
    # get the distance to the other side of the split
    nv = [0]*n
    for i in selection:
        nv[i] = v[i] + mb
    for i in complement:
        nv[i] = v[i] + ma
    return nv

def merge_trees(tree_a, tree_b):
    """
    Merge the newick trees by joining the branches with the highest serial numbers.
    The new branch connecting the trees will have the mean length of the joined branches.
    @param tree_a: a newick tree with tips marked with serial numbers
    @param tree_b: another newick tree with tips marked with serial numbers
    @return: a combined newick tree
    """
    # for each tree find the node with the highest serial number
    serial_tip_pairs = [(p.serial_number, p) for p in tree_a.gen_tips()]
    tip_a = max(serial_tip_pairs)[1]
    serial_tip_pairs = [(p.serial_number, p) for p in tree_b.gen_tips()]
    tip_b = max(serial_tip_pairs)[1]
    # reroot the trees
    tree_a.reroot(tip_a)
    tree_b.reroot(tip_b)
    # calculate the length of the new connecting branch
    neo_blen = (tree_a.root.children[0].blen + tree_b.root.children[0].blen) / 2.0
    # merge the trees
    neo_root = tree_a.root.children[0]
    neo_root.parent = None
    neo_root.blen = None
    neo_sink = tree_b.root.children[0]
    neo_sink.blen = neo_blen
    neo_root.add_child(neo_sink)
    neo_sink.set_parent(neo_root)
    # return the merged trees
    return Newick.NewickTree(neo_root)


class TreeBuilder:

    def __init__(self, initial_distance_matrix, ordered_labels, distance_matrix_splitter):
        """
        @param initial_distance_matrix: the distance matrix
        @param ordered_labels: the names of the leaves in the order they appear in the distance matrix
        @param distance_matrix_splitter: an object that can get an index selection given a distance matrix
        """
        if len(initial_distance_matrix) != len(ordered_labels):
            raise NeighborhoodJoiningError('expected a label for each distance matrix entry')
        self.initial_distance_matrix = initial_distance_matrix
        self.ordered_labels = ordered_labels
        self.splitter = distance_matrix_splitter
        self.fallback_name = 'nj'
        self.callback = None

    def set_recursion_callback(self, callback):
        """
        Set a callback function that is called every time the recursive tree building function is called.
        One application of this callback is to limit the time spent building the tree.
        The callback function takes an optional recursion depth argument.
        @param callback: None or a callable
        """
        self.callback = callback
    
    def set_fallback_name(self, fallback_name):
        """
        @param fallback_name: the name of the fallback operation to use
        """
        valid_fallback_names = ('nj', 'halving')
        if fallback_name not in valid_fallback_names:
            raise NeighborhoodJoiningError('invalid fallback name: ' + fallback_name)
        self.fallback_name = fallback_name

    def get_complexity(self):
        """
        @return: a value proportional to how long it would take to build the tree
        """
        n = len(self.ordered_labels)
        # The factor of n is due to the n splits of the tree.
        # The remaining factor is due to the time it takes to decide each split.
        return n * self.splitter.get_complexity(n)

    def build(self):
        """
        @return: the newick tree
        """
        n = len(self.ordered_labels)
        self.number_generator = SequentialNumberGenerator(n)
        self.serial_number_to_tip_set = dict((i, set([i])) for i in range(n))
        index_to_serial = range(n)
        tree = self._make_tree_helper(self.initial_distance_matrix, index_to_serial)
        # convert the serial numbers into names if possible
        for node in tree.gen_tips():
            node.add_name(self.ordered_labels[node.serial_number])
        # remove the serial numbers of all of the nodes
        for node in tree.gen_tips():
            del node.serial_number
        # convert the tree to a FelTree
        tree_string = tree.get_newick_string()
        feltree = NewickIO.parse(tree_string, FelTree.NewickTree)
        return feltree

    def on_custom_split(self, name_selection, min_split_size, max_split_size):
        """
        This is called when the custom distance matrix splitter makes a split.
        Override this callback if this information would be useful.
        @param name_selection: a set of selected names
        @param min_split_size: the size of the smaller split made by the splitter
        @param max_split_size: the size of the larger split made by the splitter
        """
        pass

    def on_nj_fallback_split(self, name_selection, min_split_size, max_split_size):
        """
        This is called when the fallback splitting method makes a split.
        This happens only immediately after the custom splitter makes a degenerate split.
        Override this callback if this information would be useful.
        @param name_selection: a set of selected names
        @param min_split_size: the size of the smaller split made by the splitter
        @param max_split_size: the size of the larger split made by the splitter
        """
        pass

    def on_halving_fallback_split(self, name_selection, min_split_size, max_split_size, halving_iteration_count):
        """
        This is called when the fallback splitting method makes a split.
        This happens only immediately after the custom splitter makes a degenerate split.
        Override this callback if this information would be useful.
        @param name_selection: a set of selected names
        @param min_split_size: the size of the smaller split made by the splitter
        @param max_split_size: the size of the larger split made by the splitter
        @param halving_iteration_count: the number of times a leaf stem length was halved
        """
        pass

    def _make_tree_helper(self, D, index_to_serial, depth=0):
        """
        Recursively build a newick tree from a distance matrix.
        @param D: a row major distance matrix
        @param index_to_serial: converts an index in D to a serial number for the tree node
        @param depth: gives the recursion depth; this is for instrumentation
        @return: a newick tree with branch lengths
        """
        # instrumentation to notify the framework that a recursive call has been made
        if self.callback:
            self.callback(depth)
        # recursively build the newick tree
        n = len(D)
        if n == 3:
            # if there are only three nodes then return a single star tree
            v, (f, g) = NeighborJoining.do_iteration(D)
            root = Newick.NewickNode()
            for i, d in enumerate(v):
                neo = Newick.NewickNode()
                neo.serial_number = index_to_serial[i]
                neo.blen = d
                root.add_child(neo)
                neo.set_parent(root)
            return Newick.NewickTree(root)
        # try to get the selection using a custom splitter
        selection = self.splitter.get_selection(D)
        complement = set(range(n)) - selection
        # if the split was insufficient then resort to either modifying the distance matrix or using neighbor joining
        fallback = False
        if min(len(selection), len(complement)) < 2:
            fallback = True
            if self.fallback_name == 'nj':
                # use an iteration of neighbor joining if this is the preferred fallback method
                v, (f, g) = NeighborJoining.do_iteration(D)
                selection = set((f, g))
                complement = set(range(n)) - selection
            elif self.fallback_name == 'halving':
                # repeatedly modify the distance matrix if this is the preferred fallback method
                halving_count = 0
                while min(len(selection), len(complement)) < 2:
                    # kill the loop if the halving count is ridiculous
                    if halving_count > 1000:
                        error_out = StringIO()
                        print >> error_out, 'the number of leaf stem halving iterations is ridiculous (%d);' % halving_count
                        print >> error_out, 'the singleton leaf stem length is %s;' % leaf_stem_length
                        print >> error_out, 'the distance matrix is:'
                        print >> error_out, MatrixUtil.m_to_string(D)
                        raise NeighborhoodJoiningError(error_out.getvalue().strip())
                    # find the index of the leaf singleton
                    halving_count += 1
                    if len(selection) == 1:
                        smaller = selection
                        larger = complement
                    elif len(complement) == 1:
                        smaller = complement
                        larger = selection
                    else:
                        error_out = StringIO()
                        print >> error_out, 'in the following distance matrix,'
                        print >> error_out, 'a split was so degenerate that it did not even leave a leaf stem to work with:'
                        print >> error_out, MatrixUtil.m_to_string(D)
                        raise NeighborhoodJoiningError(error_out.getvalue().strip())
                    v = get_crossing_distances(D, selection, complement)
                    # get the distance from the leaf singleton to the root of the rest of the tree
                    leaf_singleton_index = list(smaller)[0]
                    leaf_stem_length = v[leaf_singleton_index]
                    # if the leaf stem length is zero then repeatedly halving it will not help.
                    if not leaf_stem_length:
                        error_out = StringIO()
                        print >> error_out, 'the singleton leaf stem length is zero;'
                        print >> error_out, 'the number of leaf stem halving iterations performed was %d;' % halving_count
                        print >> error_out, 'the distance matrix is:'
                        print >> error_out, MatrixUtil.m_to_string(D)
                        raise NeighborhoodJoiningError(error_out.getvalue().strip())
                    # modify the distance matrix
                    for i in larger:
                        D[i][leaf_singleton_index] -= leaf_stem_length / 2
                        D[leaf_singleton_index][i] -= leaf_stem_length / 2
                    # get the selection and complement using the modified distance matrix
                    selection = self.splitter.get_selection(D)
                    complement = set(range(n)) - selection
        # define the new serial numbers for the selection and complement subtrees
        selection_serial = self.number_generator.get_next()
        complement_serial = self.number_generator.get_next()
        # for reporting purposes only,
        # store the subset of leaf serials defined by each new serial number
        for new_serial, indices in ((selection_serial, selection), (complement_serial, complement)):
            serials = set(index_to_serial[i] for i in indices)
            new_set = set()
            for serial in serials:
                new_set.update(self.serial_number_to_tip_set[serial])
            self.serial_number_to_tip_set[new_serial] = new_set
        # report the split
        flattened_selection = set(self.ordered_labels[serial] for serial in self.serial_number_to_tip_set[selection_serial])
        if fallback:
            if self.fallback_name == 'nj':
                self.on_nj_fallback_split(flattened_selection, len(selection), len(complement))
            elif self.fallback_name == 'halving':
                self.on_halving_fallback_split(flattened_selection, len(selection), len(complement), halving_count)
            else:
                assert False, 'internal error: invalid fallback method'
        else:
            self.on_custom_split(flattened_selection, len(selection), len(complement))
        # break the distance matrix into two distance matrices,
        # then make a tree for each one.
        A = list(sorted(selection))
        B = list(sorted(complement))
        A_distance_matrix, B_distance_matrix = split_distance_matrix(D, selection, complement)
        # define the serial numbers for the split distance matrices
        A_index_to_serial = [index_to_serial[i] for i in A] + [complement_serial]
        B_index_to_serial = [index_to_serial[i] for i in B] + [selection_serial]
        # make the next trees
        A_tree = self._make_tree_helper(A_distance_matrix, A_index_to_serial, depth+1)
        B_tree = self._make_tree_helper(B_distance_matrix, B_index_to_serial, depth+1)
        # return the merged tree
        return merge_trees(A_tree, B_tree)

def _get_validity_string(is_valid):
    if is_valid:
        return 'valid'
    else:
        return 'invalid'

class ValidatingTreeBuilder(TreeBuilder):

    def __init__(self, initial_distance_matrix, ordered_labels, distance_matrix_splitter):
        """
        @param initial_distance_matrix: the distance matrix
        @param ordered_labels: the names of the leaves in the order they appear in the distance matrix
        @param distance_matrix_splitter: an object that can get an index selection given a distance matrix
        """
        TreeBuilder.__init__(self, initial_distance_matrix, ordered_labels, distance_matrix_splitter)
        self.original_tree = None
        self.out = None
        self.mismatch_count = 0
        self.weighted_mismatch_count = 0
        self.max_mismatch_count = 0
        self.max_weighted_mismatch_count = 0

    def set_original_tree(self, tree):
        self.original_tree = tree

    def set_output_stream(self, out):
        self.out = out

    def get_weighted_mismatch_count(self):
        """
        @return: the number of incompatible partitions weighted by an n choose k measure of their depth
        """
        return self.weighted_mismatch_count

    def get_mismatch_count(self):
        """
        @return: the number of inferred tree partitions that are incompatible with the original tree
        """
        return self.mismatch_count

    def validate(self, name_selection):
        tip_selection = [tip for tip in self.original_tree.gen_tips() if tip.get_name() in name_selection]
        if self.original_tree.get_split_branch(tip_selection):
            return True
        else:
            return False

    def _on_split_helper(self, name_selection, min_split_size, max_split_size):
        """
        Do stuff that must be done on each split, and then return whether the split was valid.
        @return: True if the split was valid
        """
        # get the loss that would be incurred if the split were wrong
        n = min_split_size + max_split_size
        k = min_split_size
        mismatch_weight = Util.choose(n, k)
        # add the potential losses to the running total maximum possible loss
        if min_split_size > 1:
            self.max_mismatch_count += 1
            self.max_weighted_mismatch_count += mismatch_weight
        # if the split is actually wrong then add the losses to the running loss total
        validity = self.validate(name_selection)
        if not validity:
            self.mismatch_count += 1
            self.weighted_mismatch_count += mismatch_weight
        return validity

    def on_custom_split(self, name_selection, min_split_size, max_split_size):
        validity = self._on_split_helper(name_selection, min_split_size, max_split_size)
        validity_string = _get_validity_string(validity)
        if self.out:
            print >> self.out, '%d:%d' % (min_split_size, max_split_size), validity_string

    def on_nj_fallback_split(self, name_selection, min_split_size, max_split_size):
        validity = self._on_split_helper(name_selection, min_split_size, max_split_size)
        validity_string = _get_validity_string(validity)
        if self.out:
            print >> self.out, '%d:%d' % (min_split_size, max_split_size), validity_string, 'after falling back to neighbor joining'

    def on_halving_fallback_split(self, name_selection, min_split_size, max_split_size, halving_iterations):
        validity = self._on_split_helper(name_selection, min_split_size, max_split_size)
        validity_string = _get_validity_string(validity)
        if halving_iterations == 1:
            if self.out:
                print >> self.out, '%d:%d' % (min_split_size, max_split_size), validity_string, 'after halving the leaf stem length'
        elif halving_iterations > 1:
            if self.out:
                print >> self.out, '%d:%d' % (min_split_size, max_split_size), validity_string, 'after', halving_iterations, 'iterations of halving a leaf stem length'
        else:
            raise NeighborhoodJoiningError('expected the number of iterations to be at least one, but observed %s' % halving_iterations)


class DMSamplerError(Exception):
    """
    This exception is raised when it is apparent that the sampling will take too long, for example.
    """
    pass


class DMSampler:
    """
    Sample estimated distance matrices, rejecting degenerate distance matrices.
    Here a distance matrix is degenerate if it has entries that are estimated to
    be either zero or infinity.
    The complexity of generating the samples is also estimated,
    which is important because otherwise the rejection sampling could get stuck in a loop
    if every sample is rejected.
    """

    def __init__(self, tree, ordered_names, sequence_length):
        """
        @param tree: a phylogenetic tree object
        @param ordered_names: an ordered list of names of terminal taxa
        @param sequence_length: the length of the sequences to generate (or inf)
        """
        assert set(node.name for node in tree.gen_tips()) == set(ordered_names)
        self.tree = tree
        self.ordered_names = ordered_names
        self.sequence_length = sequence_length
        self.requested_matrix_count = 0
        self.accepted_sample_count = 0
        # reject branches of zero and infinite length
        self.zero_replacement = None
        self.inf_replacement = None
        # count the number of raw samples that had at least one zero or inf value respectively
        self.raw_zero_sample_count = 0
        self.raw_inf_sample_count = 0

    def set_zero_replacement(self, zero_replacement):
        """
        @param zero_replacement: the value replacing a zero in the distance matrix or None to reject
        """
        self.zero_replacement = zero_replacement

    def set_inf_replacement(self, inf_replacement):
        """
        @param inf_replacement: the value replacing infinity in the distance matrix or None to reject
        """
        self.inf_replacement = inf_replacement

    def _get_error_message(self):
        """
        Give a reason for failure.
        @return: an error message
        """
        # if the run has finished then report None
        if self.accepted_sample_count == self.requested_matrix_count:
            return None
        # if the run has failed when we asked for infinite sequence length then say something
        if self.sequence_length == float('inf'):
            return 'the true distance matrix was degenerate'
        # get some summary statistics of the run
        accepted = self.accepted_sample_count
        rejected_inf = self.rejected_inf_sample_count
        rejected_zero = self.rejected_zero_sample_count
        rejected = rejected_inf + rejected_zero
        total = accepted + rejected
        # otherwise make some excuse for failure
        parenthetical_remark = None
        if total:
            parenthetical_remark = '%d of %d samples accepted' % (accepted, total)
            if accepted < rejected_inf:
                parenthetical_remark += '; use shorter branch lengths or longer sequences'
            elif accepted < rejected_zero:
                parenthetical_remark += '; use longer branch lengths or longer sequences'
        error_message = 'the distance matrix sampling procedure takes too long for these settings'
        if parenthetical_remark:
            error_message += ' (%s)' % parenthetical_remark
        return error_message

    def get_rejected_sample_count(self):
        return self.rejected_inf_sample_count + self.rejected_zero_sample_count

    def get_acceptance_probability(self):
        """
        This is for progress bar stuff.
        @return: an optimistic acceptance probability for the rejection sampling
        """
        total_samples = self.accepted_sample_count + self.get_rejected_sample_count()
        if total_samples < 100:
            # if not enough samples have been taken then be optimistic
            return 1.0
        else:
            # if a reasonable number of samples have been taken then be realistic
            return self.accepted_sample_count / float(total_samples)

    def get_complexity(self):
        """
        This is for progress bar stuff.
        @return: the predicted total number of steps required, for some step granularity
        """
        # if all of the samples are rejected then the complexity is infinite
        acceptance_probability = self.get_acceptance_probability()
        if not acceptance_probability:
            return float('inf')
        # if there is some predicted probability of accepting a sample then make a guess
        n = len(self.ordered_names)
        steps_per_sample = n * n * self.sequence_length
        required_accepted_samples = self.requested_matrix_count
        samples_per_accepted_sample = 1.0 / acceptance_probability
        steps = steps_per_sample * samples_per_accepted_sample * required_accepted_samples
        return steps

    def gen_samples_or_none(self, count, max_steps):
        """
        Yield (ordered sequence list, distance matrix) pairs or None.
        The generator will stop if it sees that it cannot meet its goal in the allotted number of steps.
        The time between yielded results is bounded.
        @param count: the requested number of distance matrices or None for no bound
        @param max_steps: an upper bound on the number of steps allowed for the computation or None for no bound
        """
        # record the requested number of samples
        self.requested_matrix_count = count
        if self.sequence_length == float('inf'):
            # get the true distance matrix
            distance_matrix = self.tree.get_distance_matrix()
            # if any of the off-diagonal elements are wack then return straight away
            if self.reject_zero is None:
                if matrix_has_zero_off_diagonal(distance_matrix):
                    error_message = 'the true distance matrix has a zero off-diagonal entry'
                    raise DMSamplerError(error_message)
            if self.reject_inf:
                if matrix_has_inf_off_diagonal(distance_matrix):
                    error_message = 'the true distance matrix has an infinite off-diagonal entry'
                    raise DMSamplerError(error_message)
            # yield a bunch of copies of the true distance matrix
            for i in range(count):
                self.accepted_sample_count += 1
                yield (None, distance_matrix)
        else:
            # do some rejection sampling
            while True:
                # if we are done sampling then return
                if count is not None:
                    if self.accepted_sample_count >= count:
                        return
                # if we are taking too many computrons then bail with an error message
                if max_steps is not None:
                    if self.get_complexity() > max_steps:
                        raise DMSamplerError(self._get_error_message())
                # do the sampling
                sequence_list = JC69.sample_sequences(self.tree, self.ordered_names, self.sequence_length)
                # get the estimated distance matrix
                distance_matrix = JC69.get_ML_distance_matrix(sequence_list)
                # look for degeneracies
                if self.reject_zero and matrix_has_zero_off_diagonal(distance_matrix):
                    self.rejected_zero_sample_count += 1
                    yield None
                elif self.reject_inf and matrix_has_inf_off_diagonal(distance_matrix):
                    self.rejected_inf_sample_count += 1
                    yield None
                else:
                    self.accepted_sample_count += 1
                    yield sequence_list, distance_matrix
        

    def gen_samples(self, count, max_steps):
        """
        Yield (ordered sequence list, distance matrix) pairs.
        The generator will stop if it sees that it cannot meet its goal in the allotted number of steps.
        If no maximum number of steps is specified then the time between yielded results could be unbounded.
        @param count: the requested number of distance matrices or None for no bound
        @param max_steps: an upper bound on the number of steps allowed for the computation or None for no bound
        """
        for result in self.gen_samples_or_none(count, max_steps):
            if result is not None:
                yield result


class TestNeighborhoodJoining(unittest.TestCase):

    def test_neighbor_joining(self):
        """
        UNFINISHED
        Assert that neighbor joining is a special case of neighborhood joining.
        """
        pass


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestNeighborhoodJoining)
    unittest.TextTestRunner(verbosity=2).run(suite)

