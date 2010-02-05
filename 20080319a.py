"""Compare methods of estimating the stationary codon distribution. [FLAWED]
"""

from StringIO import StringIO
import math

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
    # define some default strings
    default_nt_string = '\n'.join(nt + ' : 1' for nt in nt_letters)
    default_aa_string = '\n'.join(aa + ' : 1' for aa in aa_letters)
    # define the form objects
    form_objects = [
            Form.MultiLine('nucleotides', 'nucleotide stationary distribution', default_nt_string),
            Form.MultiLine('aminoacids', 'amino acid stationary distribution', default_aa_string),
            Form.RadioGroup('method', 'codon stationary distribution estimation method', [
                Form.RadioItem('hb', 'equation (14) in Halpern-Bruno', True),
                Form.RadioItem('corrected', 'a corrected method')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the nucleotide distribution
    nt_to_weight = SnippetUtil.get_distribution(fs.nucleotides, 'nucleotide', nt_letters)
    # get the amino acid distribution
    aa_to_weight = SnippetUtil.get_distribution(fs.aminoacids, 'amino acid', aa_letters)
    # get distributions in convenient list form
    stationary_nt_distribution = [nt_to_weight[nt] for nt in nt_letters]
    aa_distribution = [aa_to_weight[aa] for aa in aa_letters]
    codon_distribution = []
    implied_stationary_nt_distribution = []
    if fs.corrected:
        # define the objective function
        objective_function = MyObjective(aa_distribution, stationary_nt_distribution)
        initial_guess = (0, 0, 0)
        iterations = 20
        best = scipy.optimize.nonlin.broyden2(objective_function, initial_guess, iterations)
        x, y, z = best
        best_mutation_distribution = normalized((1, math.exp(x), math.exp(y), math.exp(z)))
        # given the mutation distribution and the amino acid distribution, get the stationary distribution
        result = DirectProtein.get_nt_distribution_and_aa_energies(best_mutation_distribution, aa_distribution)
        implied_stationary_nt_distribution, result_aa_energies = result
        # get the codon distribution; kappa doesn't matter because we are only concerned with stationary distributions
        kappa = 1.0
        dpm = DirectProtein.DirectProteinRateMatrix(kappa, best_mutation_distribution, result_aa_energies)
        codon_distribution = dpm.get_stationary_distribution()
    elif fs.hb:
        # get the codon distribution
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
    # start the output text string
    out = StringIO()
    # write the codon stationary distribution
    print >> out, 'estimated codon stationary distribution:'
    for codon, p in zip(codons, codon_distribution):
        print >> out, '%s : %s' % (codon, p)
    print >> out, ''
    # write the nucleotide stationary distribution
    print >> out, 'implied nucleotide stationary distribution:'
    for nt, p in zip(nt_letters, implied_stationary_nt_distribution):
        print >> out, '%s : %s' % (nt, p)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()

def normalized(distribution):
    """
    @param distribution: a distribution
    @return: a normalized distribution
    """
    # allow the input distribution to be a generator
    P = list(distribution)
    total_weight = sum(P)
    return [p / float(total_weight) for p in P]


class MyObjective:
    """
    This is an objective function object.
    """

    def __init__(self, aa_stationary_distribution, nt_stationary_distribution):
        """
        @param aa_stationary_distribution: the amino acid stationary distribution
        @param nt_stationary_distribution: the nucleotide stationary distribution
        """
        self.aa_dist = aa_stationary_distribution
        self.nt_dist = nt_stationary_distribution
        self.guesses = []

    def __call__(self, X):
        """
        @param X: the three variables that define the mutation process nucleotide distribution.
        @return: a tuple of values that is the zero vector when the guess was right.
        """
        self.guesses.append(X)
        if len(X) != 3:
            raise ValueError('incorrect number of parameters')
        x, y, z = X
        mutation_nt_dist = normalized((1, math.exp(x), math.exp(y), math.exp(z)))
        stationary_nt_dist, aa_energies = DirectProtein.get_nt_distribution_and_aa_energies(mutation_nt_dist, self.aa_dist)
        return tuple(math.log(a/b) for a, b in zip(self.nt_dist[:3], stationary_nt_dist[:3]))
