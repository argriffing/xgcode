"""Do mixture models.
"""

# standard modules
import StringIO
import unittest

# nonstandard modules
import Codon
import MatrixUtil
import RateMatrix
import Util


class MixtureModelError(Exception):
    pass


class MixtureModel:
    """
    This class was created with a mixture of Direct Protein codon models in mind.
    """

    def __init__(self, mixture_parameters, rate_matrices):
        """
        @param mixture_parameters: an ordered list of probabilities summing to one
        @param rate_matrices: an ordered list of component rate matrix objects
        """
        if len(mixture_parameters) != len(rate_matrices):
            raise MixtureModelError('the number of mixture parameters should be the same as the number of rate matices')
        eps = .000001
        if abs(sum(mixture_parameters) - 1.0) > eps:
            raise MixtureModelError('the mixture parameters should sum to one')
        for p in mixture_parameters:
            if p < 0:
                raise MixtureModelError('no mixture parameter should be negative')
        self.mixture_parameters = mixture_parameters
        self.rate_matrices = rate_matrices

    def rescale(self, scaling_factor):
        """
        Rescale the rate matrices.
        @param scaling_factor: the factor by which to scale each component matrix
        """
        for matrix in self.rate_matrices:
            matrix.rescale(scaling_factor)

    def normalize(self):
        """
        Rescale the rate matrices so that the mixture has an expected rate of one transition per unit time.
        """
        scaling_factor = 1.0 / self.get_expected_rate()
        self.rescale(scaling_factor)
        eps = .000001
        if abs(self.get_expected_rate() - 1.0) > eps:
            raise MixtureModelError('the expected rate should be 1.0 after normalization')

    def get_expected_rate(self):
        """
        @return: the expected rate of the substitution model
        """
        expected_rates = [matrix.get_expected_rate() for matrix in self.rate_matrices]
        return Util.dot_product(self.mixture_parameters, expected_rates)

    def get_stationary_distribution(self):
        stationary_distributions = [matrix.get_stationary_distribution() for matrix in self.rate_matrices]
        return [Util.dot_product(proportions, self.mixture_parameters) for proportions in zip(*stationary_distributions)]

    def simulate_states(self, tree):
        """
        Decorate a tree by adding simulated states to every node.
        @param tree: a tree with branch lengths
        """
        weight_matrix_pairs = zip(self.mixture_parameters, self.rate_matrices)
        matrix = Util.weighted_choice(weight_matrix_pairs)
        matrix.simulate_states(tree)

    def simulate_ancestral_states(self, tree):
        """
        Decorate a tree by adding simulated states conditional on its leaf states.
        @param tree: a tree with branch lengths and leaf states
        """
        conditional_mixture_parameters = self.get_membership(tree)
        weight_matrix_pairs = zip(conditional_mixture_parameters, self.rate_matrices)
        matrix = Util.weighted_choice(weight_matrix_pairs)
        matrix.simulate_ancestral_states(tree)

    def get_likelihood(self, tree):
        """
        Calculate the likelihood.
        @param tree: a tree with branch lengths and leaf states
        """
        likelihoods = [matrix.get_likelihood(tree) for matrix in self.rate_matrices]
        return Util.dot_product(self.mixture_parameters, likelihoods)

    def get_membership(self, tree):
        """
        Calculate the maximum likelihood membership of the site.
        @param tree: a tree with branch lengths and leaf states
        @return: an ordered list of membership proportions
        """
        likelihoods = [matrix.get_likelihood(tree) for matrix in self.rate_matrices]
        total = Util.dot_product(self.mixture_parameters, likelihoods)
        return [(p * likelihood) / total for p, likelihood in zip(self.mixture_parameters, likelihoods)]


class TestSubModel(unittest.TestCase):

    def test_foo(self):
        pass


def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestRateMatrix)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()


