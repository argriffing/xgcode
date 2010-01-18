
import unittest
import StringIO
import math
import random

import Util
import Newick
import Fasta
import RateMatrix
import MatrixUtil

class SimulationError(Exception):
    pass

def get_log_likelihood(tree, alignment, substitution_model):
    """
    @param tree: a newick tree with branch lengths
    @param alignment: a Fasta Alignment object with headers that match the tree tip names
    @param substitution_model: a way to get a likelihood from a tree given its leaf states
    @return: the log likelihood or None if there is no likelihood
    """
    # Get the number of times each different column appears.
    column_multiset = alignment.get_column_multiset()
    # Get the likelihood.
    log_likelihood = 0
    for col, count in column_multiset.items():
        name_to_state = dict(zip(alignment.headers, col))
        # Augment each tip with its corresponding letter.
        for tip in tree.gen_tips():
            tip.state = name_to_state[tip.name]
        # Add the contribution of the current column pattern to the total log likelihood.
        column_likelihood = substitution_model.get_likelihood(tree)
        if not column_likelihood:
            return None
        log_likelihood += math.log(column_likelihood) * count
    return log_likelihood

def simulate_alignment(tree, substitution_model, ncolumns, seed=None):
    """
    @param tree: a newick tree with branch lengths
    @param substitution_model: a way to simulate states on a tree
    @param ncolumns: the number of columns to simulate
    @param seed: a random number seed
    @return: a Fasta Alignment object of the simulated sequences
    """
    # Check the input.
    for node in tree.gen_non_root_nodes():
        if node.get_branch_length() is None or node.get_branch_length() <= 0:
            raise SimulationError('all branch lengths should be positive')
    tip_names = [node.name for node in tree.gen_tips()]
    for name in tip_names:
        if not name:
            raise SimulationError('each leaf should have a name')
    if len(tip_names) != len(set(tip_names)):
        raise SimulationError('each leaf should have a unique name')
    # Save the rng state if we are using a seed.
    if seed is not None:
        old_rng_state = random.getstate()
    # Seed the rng if we are using a seed.
    if seed is not None:
        random.seed(seed)
    # Simulate the states on the tree.
    simulated_sequences = dict((node.name, []) for node in tree.gen_tips())
    for column_index in range(ncolumns):
        substitution_model.simulate_states(tree)
        for node in tree.gen_tips():
            simulated_sequences[node.name].append(node.state)
    # Restore the rng state if we are using a seed
    if seed is not None:
        random.setstate(old_rng_state)
    # Create an alignment object from the simulated sequences.
    sio = StringIO.StringIO()
    for header, sequence in simulated_sequences.items():
        print >> sio, '>' + header
        print >> sio, ''.join(sequence)
    fasta_string = sio.getvalue()
    return Fasta.Alignment(StringIO.StringIO(fasta_string))

def simulate_ancestral_alignment(tree, alignment, substitution_model):
    """
    @param tree: a newick tree with branch lengths
    @param alignment: a Fasta Alignment object with headers that match the tree tip names
    @param substitution_model: a way to simulate ancestral states from a tree given its leaf states
    @return: a Fasta Alignment object of the simulated ancestral sequences
    """
    for node in tree.gen_non_root_nodes():
        if node.get_branch_length() is None or node.get_branch_length() <= 0:
            raise SimulationError('all branch lengths should be positive')
    for node in tree.gen_internal_nodes():
        if not node.name:
            raise SimulationError('all internal nodes should be named')
    simulated_ancestors = dict((node.name, []) for node in tree.gen_internal_nodes())
    for col in alignment.columns:
        name_to_letter = dict(zip(alignment.headers, col))
        # Augment each tip with its corresponding letter.
        for tip in tree.gen_tips():
            tip.state = name_to_letter[tip.name]
        # Do the simulation.
        substitution_model.simulate_ancestral_states(tree)
        name_state_pairs = [(node.name, node.state) for node in tree.gen_internal_nodes_preorder()]
        # Add this simulated column.
        for name, state in name_state_pairs:
            simulated_ancestors[name].append(state)
    # Create an alignment object from the simulated sequences.
    sio = StringIO.StringIO()
    print >> sio, alignment.to_fasta_string()
    for header, sequence in simulated_ancestors.items():
        print >> sio, '>' + header
        print >> sio, ''.join(sequence)
    fasta_string = sio.getvalue()
    return Fasta.Alignment(StringIO.StringIO(fasta_string))


class TestPhyLikelihood(unittest.TestCase):

    def test_likelihood(self):
        # Parse the example tree.
        tree_string = Newick.brown_example_tree
        tree = Newick.parse(tree_string, Newick.NewickTree)
        tree.assert_valid()
        # Get header and sequence pairs.
        alignment = Fasta.Alignment(StringIO.StringIO(Fasta.brown_example_alignment))
        # Get the Jukes-Cantor rate matrix object.
        dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
        ordered_states = list('ACGT')
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, ordered_states, ordered_states)
        rate_matrix_object = RateMatrix.RateMatrix(row_major_rate_matrix, ordered_states)
        # Calculate the log likelihood.
        log_likelihood = get_log_likelihood(tree, alignment, rate_matrix_object)
        self.assertAlmostEqual(log_likelihood, -4146.26547208)

    def test_simulation(self):
        tree_string = '(((Human:0.1, Chimpanzee:0.2)to-chimp:0.8, Gorilla:0.3)to-gorilla:0.7, Orangutan:0.4, Gibbon:0.5)all;'
        # Parse the example tree.
        tree = Newick.parse(tree_string, Newick.NewickTree)
        tree.assert_valid()
        # Get header and sequence pairs.
        alignment = Fasta.Alignment(StringIO.StringIO(Fasta.brown_example_alignment))
        # Get the Jukes-Cantor rate matrix object.
        dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix()
        ordered_states = list('ACGT')
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, ordered_states, ordered_states)
        rate_matrix_object = RateMatrix.RateMatrix(row_major_rate_matrix, ordered_states)
        # Simulate ancestral states.
        simulated_alignment = simulate_ancestral_alignment(tree, alignment, rate_matrix_object)


def main():
    print 'try --test'


if __name__ == '__main__':
    from optparse import OptionParser
    import profile
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPhyLikelihood)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

