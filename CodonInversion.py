"""This module demonstrates a problem I had with broyden2.
"""

import unittest
import math
import scipy.optimize
import logging


# nucleotide letters
nt_letters = ('A', 'C', 'G', 'T')

# amino acid letters
aa_letters = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')

# standard codons that are not stop codons
codons = (
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG',
        'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC',
        'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA',
        'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT',
        'GTA', 'GTC', 'GTG', 'GTT', 'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC',
        'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'
        )

# define the standard genetic code neglecting stop codons
codon_to_aa_letter = {
        'ACC': 'T', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 'ATC': 'I', 'AAC': 'N',
        'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'CTC': 'L', 'AGC': 'S', 'ACA': 'T',
        'CTT': 'L', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L',
        'ACT': 'T', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'AGT': 'S', 'CAG': 'Q',
        'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'TGT': 'C', 'CGA': 'R',
        'CCA': 'P', 'TCT': 'S', 'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C',
        'GGG': 'G', 'GGA': 'G', 'TGG': 'W', 'GGC': 'G', 'TAC': 'Y', 'TTC': 'F',
        'TCG': 'S', 'TTA': 'L', 'AGA': 'R', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E',
        'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 'GTC': 'V', 'GCG': 'A',
        'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GCT': 'A', 'TTG': 'L', 'CGT': 'R',
        'CGC': 'R'
        }

# for sanity checking only
eps = .000000001

def almost_equals(a, b):
    return abs(a-b) < eps

def normalized(distribution):
    """
    @param distribution: a sequence of weights
    @return: a normalized distribution
    """
    # allow the input distribution to be a generator
    P = list(distribution)
    total_weight = sum(P)
    return [p / float(total_weight) for p in P]

def product(numbers):
    x = 1
    for number in numbers:
        x *= number
    return x


class Objective:
    """
    An objective function that stores guesses and responses.
    """

    def __init__(self):
        self.guesses = []
        self.responses = []

    def __call__(self, X):
        """
        This is called by the optimization framework.
        """
        self.guesses.append(X)
        try:
            response = self.evaluate(X)
        except Exception, e:
            self.responses.append('???')
            raise e
        self.responses.append(response)
        return response

    def log_history(self):
        for guess, response in zip(self.guesses, self.responses):
            logging.debug('%s -> %s', str(guess), str(response))

    def evaluate(self, X):
        raise NotImplementedError('override this member function')


class CodonObjective(Objective):
    """
    This is an objective function.
    """

    def __init__(self, aa_stationary_distribution, nt_stationary_distribution):
        """
        @param aa_stationary_distribution: the amino acid stationary distribution
        @param nt_stationary_distribution: the nucleotide stationary distribution
        """
        Objective.__init__(self)
        self.aa_dist = aa_stationary_distribution
        self.nt_dist = nt_stationary_distribution

    def evaluate(self, X):
        """
        @param X: the three variables that define the mutation process nucleotide distribution.
        @return: a tuple of values that is the zero vector when the guess was right.
        """
        if len(X) != 3:
            raise ValueError('incorrect number of parameters')
        x, y, z = X
        mutation_nt_dist = normalized((1, math.exp(x), math.exp(y), math.exp(z)))
        if not almost_equals(sum(mutation_nt_dist), 1.0):
            raise ValueError('detected possibly invalid input to the objective function: ' + str(X))
        stationary_nt_dist, aa_energies = get_nt_distribution_and_aa_energies(mutation_nt_dist, self.aa_dist)
        return tuple(1 + math.log(a/b)**2 for a, b in zip(self.nt_dist[:3], stationary_nt_dist[:3]))


class SimpleObjective(Objective):
    """
    This is an objective function.
    """

    def evaluate(self, X):
        return tuple(x*x - 1 for x in X)
        #for i, x in enumerate(X):
            #difference = x - i
            #result = difference * (1 + difference)
            #results.append(result)
        #return results


def get_nt_distribution_and_aa_energies(mutation_distribution, aa_distribution):
    """
    Use this function to guess the mutation distribution.
    If the output of this function is near the observed nucleotide distribution,
    then the input mutation distribution was almost correct.
    @param mutation_distribution: an ordered list of nucleotide frequencies defined by the mutation process
    @param aa_distribution: an ordered list of amino acid frequencies defined by both the mutation and selection process
    @return: the observed nucleotide distribution conditioned on the input distributions, and the amino acid energies
    """
    # do some error checking
    if len(mutation_distribution) != 4:
        raise ValueError('expected four nucleotides')
    if len(aa_distribution) != 20:
        raise ValueError('expected twenty amino acids')
    if not almost_equals(sum(mutation_distribution), 1.0):
        raise ValueError('the mutation distribution does not sum to 1.0')
    if not almost_equals(sum(aa_distribution), 1.0):
        raise ValueError('the amino acid distribution does not sum to 1.0')
    for value in mutation_distribution:
        if almost_equals(value, 0):
            raise ValueError('each nucleotide should have a positive weight')
    for value in aa_distribution:
        if almost_equals(value, 0):
            raise ValueError('each amino acid should have a positive weight')
    # first get codon weights that depend only on the mutation distribution
    # and get the amino acid weights that depend only on the mutation distribution
    nt_to_weight = dict(zip(nt_letters, mutation_distribution))
    aa_to_weight = dict(zip(aa_letters, [0] * 20))
    codon_to_weight = {}
    for codon in codons:
        aa = codon_to_aa_letter[codon]
        weight = product(nt_to_weight[nt] for nt in codon)
        codon_to_weight[codon] = weight
        aa_to_weight[aa] += weight
    # rescale codon and amino acid weights to sum to one
    total_weight = sum(aa_to_weight.values())
    for codon in codons:
        codon_to_weight[codon] /= total_weight
    for aa in aa_letters:
        aa_to_weight[aa] /= total_weight
    # now find the amino acid exponentiated negative energies that scale the codons to the correct stationary distribution
    aa_to_exp_neg_energy = {}
    for aa, target_proportion in zip(aa_letters, aa_distribution):
        aa_to_exp_neg_energy[aa] = target_proportion / aa_to_weight[aa]
    # now recalculate the codon weights to match the target aa distribution
    for codon in codons:
        aa = codon_to_aa_letter[codon]
        codon_to_weight[codon] *= aa_to_exp_neg_energy[aa]
    if not almost_equals(sum(codon_to_weight.values()), 1.0):
        raise RuntimeError('the final codon weights do not sum to 1.0')
    # get the final nucleotide weights
    nt_to_final_weight = dict(zip(nt_letters, [0]*4))
    for codon, weight in codon_to_weight.items():
        for nt in codon:
            nt_to_final_weight[nt] += weight
    total_nt_weight = float(sum(nt_to_final_weight.values()))
    for nt in nt_letters:
        nt_to_final_weight[nt] /= total_nt_weight
    if not almost_equals(sum(nt_to_final_weight.values()), 1.0):
        raise RuntimeError('the final nucleotide weights do not sum to 1.0')
    # get the final nucleotide list
    final_nucleotide_list = [nt_to_final_weight[nt] for nt in nt_letters]
    # get the final amino acid list
    final_amino_acid_list = [-math.log(aa_to_exp_neg_energy[aa]) for aa in aa_letters]
    # ok center the final amino acid list
    mean_energy = sum(final_amino_acid_list) / float(len(final_amino_acid_list))
    final_amino_acid_list = [energy - mean_energy for energy in final_amino_acid_list]
    # return the lists
    return final_nucleotide_list, final_amino_acid_list


class TestCodonInversion(unittest.TestCase):

    def _test_solver_helper(self, iterations):
        # create some uniform stationary distributions
        aa_distribution = [0.05] * 20
        observed_nt_stationary_distribution = [0.25] * 4
        # define the objective function
        objective_function = CodonObjective(aa_distribution, observed_nt_stationary_distribution)
        # solve f(x, y, z) = (0, 0, 0)
        initial_guess = (0, 0, 0)
        try:
            best = scipy.optimize.nonlin.broyden2(objective_function, initial_guess, iterations)
            objective_function.log_history()
        except Exception, e:
            error = RuntimeError('failed on evaluation %d of %d: %s' % (len(objective_function.guesses), iterations, str(e)))
            objective_function.log_history()
            raise error
        x, y, z = best
        best_mutation_distribution = normalized((1, math.exp(x), math.exp(y), math.exp(z)))

    def test_solver_20(self):
        self._test_solver_helper(20)

    def test_solver_40(self):
        self._test_solver_helper(40)

    def test_solver_simple(self):
        f = SimpleObjective()
        iterations = 40
        try:
            scipy.optimize.nonlin.broyden2(f, (0, 0, 0), iterations)
            f.log_history()
        except Exception, e:
            error = RuntimeError('failed on evaluation %d of %d: %s' % (len(f.guesses), iterations, str(e)))
            f.log_history()
            raise error



def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        if options.verbose:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.WARNING)
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCodonInversion)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()
