"""Infer the mutational stationary distribution from the nucleotide and amino acid stationary distributions. [FLAWED]

Given a nucleotide stationary distribution and an amino acid stationary distribution,
infer the mutation and selection components.
The mutation component is the stationary distribution of the nucleotide mutation process.
"""

# The selection component is centered amino acid energy vector upon which selection operates.

from StringIO import StringIO
import math
import random

import scipy
import scipy.optimize

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import Codon
import DirectProtein
import Form
from Codon import g_sorted_nt_letters as nt_letters
from Codon import g_sorted_aa_letters as aa_letters
from Codon import g_sorted_non_stop_codons as codons

def get_form():
    """
    @return: the body of a form
    """
    # define the default distributions
    default_nt_string = '\n'.join(nt + ' : 1' for nt in sorted(Codon.g_nt_letters))
    default_aa_string = '\n'.join(aa + ' : 1' for aa in sorted(Codon.g_aa_letters))
    # define the form objects
    form_objects = [
            Form.MultiLine('nucleotides', 'nucleotide stationary distribution', default_nt_string),
            Form.MultiLine('aminoacids', 'amino acid stationary distribution', default_aa_string)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the nucleotide distribution
    nt_to_probability = SnippetUtil.get_distribution(fs.nucleotides, 'nucleotide', nt_letters)
    # get the amino acid distribution
    aa_to_probability = SnippetUtil.get_distribution(fs.aminoacids, 'amino acid', aa_letters)
    # convert the dictionaries to lists
    observed_nt_stationary_distribution = [nt_to_probability[nt] for nt in nt_letters]
    aa_distribution = [aa_to_probability[aa] for aa in aa_letters]
    # define the objective function
    objective_function = MyCodonObjective(aa_distribution, observed_nt_stationary_distribution)
    initial_stationary_guess = halpern_bruno_nt_estimate(nt_to_probability, aa_to_probability)
    A, C, G, T = initial_stationary_guess
    initial_guess = (math.log(C/A), math.log(G/A), math.log(T/A))
    iterations = 20
    try:
        best = scipy.optimize.nonlin.broyden2(objective_function, initial_guess, iterations)
    except Exception, e:
        debugging_information = objective_function.get_history()
        raise HandlingError(str(e) + '\n' + debugging_information)
    x, y, z = best
    best_mutation_distribution = normalized((1, math.exp(x), math.exp(y), math.exp(z)))
    # given the mutation distribution and the amino acid distribution, get the stationary distribution
    result = DirectProtein.get_nt_distribution_and_aa_energies(best_mutation_distribution, aa_distribution)
    result_stationary_nt_dist, result_aa_energies = result
    # make a results string
    out = StringIO()
    # write the stationary nucleotide distribution of the mutation process
    print >> out, 'mutation nucleotide stationary distribution:'
    for nt, probability in zip(nt_letters, best_mutation_distribution):
        print >> out, '%s : %s' % (nt, probability)
    # write the centered amino acid energies
    print >> out, ''
    print >> out, 'amino acid energies:'
    for aa, energy in zip(aa_letters, result_aa_energies):
        print >> out, '%s : %s' % (aa, energy)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

# for sanity checking only
eps = .000000001

def almost_equals(a, b):
    return abs(a-b) < eps

def normalized(distribution):
    """
    @param distribution: a distribution
    @return: a normalized distribution
    """
    # allow the input distribution to be a generator
    P = list(distribution)
    total_weight = sum(P)
    return [p / float(total_weight) for p in P]

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

    def get_history(self):
        out = StringIO()
        for guess, response in zip(self.guesses, self.responses):
            print >> out, str(guess), '->', str(response)
        return out.getvalue()

    def evaluate(self, X):
        raise NotImplementedError('override this member function')


class CodonObjective(Objective):

    def __init__(self):
        Objective.__init__(self)
        self.energies = []

    def __call__(self, X):
        """
        This is called by the optimization framework.
        """
        self.guesses.append(X)
        try:
            response, energies = self.evaluate(X)
        except Exception, e:
            self.responses.append('???')
            self.energies.append('???')
            raise e
        self.responses.append(response)
        self.energies.append(energies)
        energy_penalty = sum(abs(energy) for energy in energies)
        #return tuple(x + energy_penalty for x in response)
        return response

    def get_history(self):
        out = StringIO()
        for guess, response, energies in zip(self.guesses, self.responses, self.energies):
            print >> out, str(guess), '->', str(response), ':', str(energies)
        return out.getvalue()

    def evaluate(self, X):
        raise NotImplementedError('override this member function')


class MyCodonObjective(CodonObjective):
    """
    This is an objective function.
    """

    def __init__(self, aa_stationary_distribution, nt_stationary_distribution):
        """
        @param aa_stationary_distribution: the amino acid stationary distribution
        @param nt_stationary_distribution: the nucleotide stationary distribution
        """
        CodonObjective.__init__(self)
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
        stationary_nt_dist, aa_energies = DirectProtein.get_nt_distribution_and_aa_energies(mutation_nt_dist, self.aa_dist)
        evaluation = tuple(math.log(a/b) for a, b in zip(self.nt_dist[:3], stationary_nt_dist[:3]))
        return evaluation, aa_energies


def halpern_bruno_nt_estimate(nt_to_weight, aa_to_weight):
    """
    @param nt_to_weight: a dictionary specifying the nucleotide stationary distribution
    @param aa_to_weight: a dictionary specifying the amino acid stationary distribution
    @return: a vector of estimated nucleotide proportions
    """
    unnormalized_codon_distribution = []
    for codon in codons:
        aa = Codon.g_codon_to_aa_letter[codon]
        sibling_codons = Codon.g_aa_letter_to_codons[aa]
        codon_aa_weight = aa_to_weight[aa]
        codon_nt_weight = Util.product(nt_to_weight[nt] for nt in codon)
        sibling_nt_weight_sum = sum(Util.product(nt_to_weight[nt] for nt in sibling) for sibling in sibling_codons)
        weight = (codon_aa_weight * codon_nt_weight) / sibling_nt_weight_sum
        unnormalized_codon_distribution.append(weight)
    codon_distribution = normalized(unnormalized_codon_distribution)
    nt_to_weight = dict(zip(nt_letters, [0]*4))
    for codon, p in zip(codons, codon_distribution):
        for nt in codon:
            nt_to_weight[nt] += p
    implied_stationary_nt_distribution = normalized(nt_to_weight[nt] for nt in nt_letters)
    return implied_stationary_nt_distribution 
