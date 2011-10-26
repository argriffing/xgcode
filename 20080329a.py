"""Sample a nucleotide alignment given a tree and a HKY mixture.

The mixture is scaled so that the branch lengths in the newick tree
are the expected number of substitutions on the branch.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import RateMatrix
import Newick
import PhyLikelihood
import MatrixUtil
import SubModel
import Nexus
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.Integer('ncols', 'sample this many nucleotide columns',
                100, low=1, high=1000),
            Form.MultiLine('frequency_a', 'first component frequencies',
                get_frequency_string(0)),
            Form.Float('kappa_a', 'first component kappa',
                get_kappa(0), low_inclusive=0),
            Form.Float('weight_a', 'first component weight',
                get_weight(0), low_inclusive=0),
            Form.MultiLine('frequency_b', 'second component frequencies',
                get_frequency_string(1)),
            Form.Float('kappa_b', 'second component kappa',
                get_kappa(1), low_inclusive=0),
            Form.Float('weight_b', 'second component weight',
                get_weight(1), low_inclusive=0),
            Form.RadioGroup('fmt', 'output format options', [
                Form.RadioItem('fasta', 'fasta'),
                Form.RadioItem('nex', 'nexus', True)]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Alignment('out.%s', interpolants=['fmt'])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # parse the tree
    try:
        tree = Newick.parse(fs.tree, Newick.NewickTree)
        tree.assert_valid()
    except Newick.NewickSyntaxError, e:
        raise HandlingError(str(e))
    # get the mixture weights
    mixture_weights = [fs.weight_a, fs.weight_b]
    # get the kappa values
    kappa_values = [fs.kappa_a, fs.kappa_b]
    # get the nucleotide distributions
    frequency_strings = (fs.frequency_a, fs.frequency_b)
    nucleotide_distributions = []
    for nt_string in frequency_strings:
        d = SnippetUtil.get_distribution(nt_string, 'nucleotide', list('ACGT'))
        nucleotide_distributions.append(d)
    # create the nucleotide HKY rate matrix objects
    rate_matrix_objects = []
    for nt_distribution, kappa in zip(nucleotide_distributions, kappa_values):
        rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(
                nt_distribution, kappa)
        rate_matrix_objects.append(rate_matrix_object)
    # create the mixture proportions
    weight_sum = sum(mixture_weights)
    mixture_proportions = [weight / weight_sum for weight in mixture_weights]
    # create the mixture model
    mixture_model = SubModel.MixtureModel(
            mixture_proportions, rate_matrix_objects)
    # normalize the mixture model
    mixture_model.normalize()
    # simulate the alignment
    try:
        alignment = PhyLikelihood.simulate_alignment(
                tree, mixture_model, fs.ncols)
    except PhyLikelihood.SimulationError, e:
        raise HandlingError(e)
    # get the output string
    output_string = ''
    if fs.fasta:
        # the output is the alignment
        arr = []
        for node in tree.gen_tips():
            arr.append(alignment.get_fasta_sequence(node.name))
        alignment_string = '\n'.join(arr)
        output_string = alignment_string
    elif fs.nex:
        # the output is the alignment and the tree
        nexus = Nexus.Nexus()
        nexus.tree = tree
        nexus.alignment = alignment
        for i in range(2):
            arr = []
            arr.append('weight: %s' % mixture_weights[i])
            arr.append('kappa: %s' % kappa_values[i])
            nexus.add_comment('category %d: %s' % (i+1, ', '.join(arr)))
        output_string = str(nexus)
    # define the filename
    if fs.fasta:
        filename_extension = 'fasta'
    elif fs.nex:
        filename_extension = 'nex'
    filename = 'sample.' + fs.fmt
    # send the response
    contentdisposition = "%s; filename=%s" % (fs.contentdisposition, filename)
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', contentdisposition)]
    return response_headers, output_string

def get_kappa(index):
    return [2, 2][index]

def get_weight(index):
    return [3, 1][index]

def gen_frequency_lines(index):
    category_to_frequencies = [
            [1, 1, 1, 1],
            [1, 4, 4, 1]]
    frequencies = category_to_frequencies[index]
    for nt, frequency in zip('ACGT', frequencies):
        yield '%s : %s' % (nt, frequency)

def get_single_line_frequency_string(index):
    return ', '.join(gen_frequency_lines(index))

def get_frequency_string(index):
    return '\n'.join(gen_frequency_lines(index))
