"""Get the maximum likelihood rate for each column of a nucleotide alignment.
"""

from StringIO import StringIO

import scipy.optimize

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import Fasta
import RateMatrix
import PhyLikelihood
import Stockholm
import Discretizer
import Form
import FormOut

def get_stockholm_string(tree, alignment, mle_rates):
    """
    @param alignment: a nucleotide alignment
    @param mle_rates: a mle rate for each column
    @return: a string that represents a file in stockholm format
    """
    # define the bins
    max_categories = 10
    discretizer = Discretizer.Discretizer(mle_rates, max_categories)
    boundaries = discretizer.get_boundaries()
    # determine to which bin each rate belongs
    categorized_rates = []
    for rate in mle_rates:
        for i, (low, high) in enumerate(boundaries):
            if low <= rate <= high:
                categorized_rates.append(i)
                continue
    # define the annotation string
    annotation_name = 'MLE_rate'
    annotation_value = ''.join(str(x) for x in categorized_rates)
    # define the stockholm object
    stockholm = Stockholm.Stockholm()
    stockholm.tree = tree
    stockholm.alignment = alignment
    # add comments and annotations to the stockholm object
    stockholm.add_comment('MLE_rate ranges:')
    for i, b in enumerate(boundaries):
        stockholm.add_comment(str(i) + ': ' + str(b))
    stockholm.add_column_annotation(annotation_name, annotation_value)
    return str(stockholm)

class Objective:
    """
    This is a function object that evaluates the likelihood of various rates.
    """

    def __init__(self, tree, rate_matrix):
        """
        @param tree: a phylogenetic tree with branch lengths and annotated tips
        @param rate_matrix: a rate matrix used to calculate column likelihoods
        """
        self.tree = tree
        self.rate_matrix = rate_matrix

    def __call__(self, rate):
        """
        Return the negative likelihood of a column.
        The negative likelihood is computed using
        the tree, matrix, and rate.
        @param rate: the rate of the rate matrix
        @return: the negative likelihood of the column
        """
        if not rate:
            inf = float('inf')
            neginf = float('-inf')
            states = [tip.state for tip in self.tree.gen_tips()]
            if len(set(states)) == 1:
                likelihood = 1
            else:
                likelihood = 0
        else:
            self.rate_matrix.set_rate(rate)
            likelihood = RateMatrix.get_likelihood(self.tree, self.rate_matrix)
        return -likelihood


def get_mle_rates(tree, alignment, rate_matrix):
    """
    @param tree: a tree with branch lengths
    @param alignment: a nucleotide alignment
    @param rate_matrix: a nucleotide rate matrix object
    @return: a list giving the maximum likelihood rate for each column
    """
    # define the objective function
    objective_function = Objective(tree, rate_matrix)
    # create the cache so each unique column is evaluated only once
    column_to_mle_rate = {}
    # look for maximum likelihood rates
    mle_rates = []
    for column in alignment.columns:
        column_tuple = tuple(column)
        if column_tuple in column_to_mle_rate:
            mle_rate = column_to_mle_rate[column_tuple]
        else:
            if len(set(column)) == 1:
                # If the column is homogeneous
                # then we know that the mle rate is zero.
                mle_rate = 0
            else:
                # redecorate the tree with nucleotide states at the tips
                name_to_state = dict(zip(alignment.headers, column))
                for tip in tree.gen_tips():
                    tip.state = name_to_state[tip.name]
                # Get the maximum likelihood rate for the column
                # using a golden section search.
                # The bracket is just a suggestion.
                bracket = (0.51, 2.01)
                mle_rate = scipy.optimize.golden(
                        objective_function, brack=bracket)
            column_to_mle_rate[column_tuple] = mle_rate
        mle_rates.append(mle_rate)
    return mle_rates

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = Newick.brown_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the default alignment string
    default_alignment_string = Fasta.brown_example_alignment.strip()
    # define the default nucleotide distribution weights
    default_nt_lines = [
            'A : 1',
            'C : 1',
            'G : 1',
            'T : 1']
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('alignment', 'nucleotide alignment',
                default_alignment_string),
            Form.MultiLine('weights', 'nucleotide weights',
                '\n'.join(default_nt_lines)),
            Form.Float('kappa', 'transition transversion rate ratio',
                2, low_inclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.Stockholm()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the nucleotide distribution
    distribution = SnippetUtil.get_distribution(
            fs.weights, 'nucleotide', list('ACGT'))
    # get the nucleotide alignment
    try:
        alignment = Fasta.Alignment(StringIO(fs.alignment))
        alignment.force_nucleotide()
    except Fasta.AlignmentError, e:
        raise HandlingError(e)
    # get the rate matrix defined by the nucleotide distribution and kappa
    row_major_rate_matrix = RateMatrix.get_unscaled_hky85_rate_matrix(
            distribution, fs.kappa).get_row_major_rate_matrix()
    rate_matrix = RateMatrix.FastRateMatrix(
            row_major_rate_matrix, list('ACGT'))
    rate_matrix.normalize()
    # get the mle rates
    mle_rates = get_mle_rates(tree, alignment, rate_matrix)
    # get the response
    stockholm_string = get_stockholm_string(tree, alignment, mle_rates)
    # return the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, stockholm_string

def main():
    # create the alignment object
    print 'creating the alignment...'
    alignment_string = Fasta.brown_example_alignment.strip()
    alignment = Fasta.Alignment(StringIO(alignment_string))
    # create a tree object
    print 'creating the tree...'
    tree_string = Newick.brown_example_tree
    tree = Newick.parse(tree_string, Newick.NewickTree)
    # create a rate matrix object
    print 'creating the rate matrix object...'
    distribution = {'A': .25, 'C': .25, 'G': .25, 'T': .25}
    kappa = 2.0
    row_major_rate_matrix = RateMatrix.get_unscaled_hky85_rate_matrix(
            distribution, kappa).get_row_major_rate_matrix()
    rate_matrix = RateMatrix.FastRateMatrix(
            row_major_rate_matrix, list('ACGT'))
    rate_matrix.normalize()
    # get the mle_rates
    print 'getting the mle rates...'
    mle_rates = get_mle_rates(tree, alignment, rate_matrix)
    print 'mle rates:'
    print mle_rates
    print 'stockholm string:'
    print get_stockholm_string(tree, alignment, mle_rates)

if __name__ == '__main__':
    #import profile
    #profile.run('main()')
    main()


