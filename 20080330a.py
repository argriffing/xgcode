"""Calculate the likelihood given nexus data and a mixture of two HKY models.

The nexus data should have a tree and an alignment.
The mixture is scaled so that the branch lengths in the newick tree are the expected number of substitutions on the branch.
"""

import math
from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import RateMatrix
import Monospace
import HeatMap
import Util
import Newick
import SubModel
import Nexus
import Form

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('nexus', 'nexus data', Nexus.nexus_sample_string.strip()),
            Form.MultiLine('frequency_a', 'first component nucleotide frequencies', get_frequency_string(0)),
            Form.Float('kappa_a', 'first component kappa', get_kappa(0), low_inclusive=0),
            Form.Float('weight_a', 'first component weight', get_weight(0), low_inclusive=0),
            Form.MultiLine('frequency_b', 'second component nucleotide frequencies', get_frequency_string(1)),
            Form.Float('kappa_b', 'second component kappa', get_kappa(1), low_inclusive=0),
            Form.Float('weight_b', 'second component weight', get_weight(1), low_inclusive=0)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the nexus data
    nexus = Nexus.Nexus()
    try:
        nexus.load(StringIO(fs.nexus))
    except Nexus.NexusError, e:
        raise HandlingError(e)
    # get the mixture weights
    mixture_weights = [fs.weight_a, fs.weight_b]
    # get the kappa values
    kappa_values = [fs.kappa_a, fs.kappa_b]
    # get the nucleotide distributions
    nucleotide_distributions = []
    for nt_string in (fs.frequency_a, fs.frequency_b):
        distribution = SnippetUtil.get_distribution(nt_string, 'nucleotide', list('ACGT'))
        nucleotide_distributions.append(distribution)
    # create the nucleotide HKY rate matrix objects
    rate_matrix_objects = []
    for nt_distribution, kappa in zip(nucleotide_distributions, kappa_values):
        rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(nt_distribution, kappa)
        rate_matrix_objects.append(rate_matrix_object)
    # create the mixture proportions
    weight_sum = sum(mixture_weights)
    mixture_proportions = [weight / weight_sum for weight in mixture_weights]
    # create the mixture model
    mixture_model = SubModel.MixtureModel(mixture_proportions, rate_matrix_objects)
    # normalize the mixture model
    mixture_model.normalize()
    # create the results
    html_string = do_analysis(mixture_model, nexus.alignment, nexus.tree)
    # write the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, html_string

def do_analysis_helper(labels, element_lists, width):
    """
    Chop up the rows of data.
    Yield lines of text to be displayed in an html pre tag.
    @param labels: row labels to be left justified
    @param element_lists: the elements of each data row where each element is a letter or a span
    @param width: the number of elements allowed per page row
    """
    if len(set(len(element_list) for element_list in element_lists)) != 1:
        raise ValueError('each element list should have the same nonzero number of elements')
    label_width = max(len(label) for label in labels) + 1
    chopped_element_lists = [list(Util.chopped(element_list, width)) for element_list in element_lists]
    page_rows = zip(*chopped_element_lists)
    for i, page_row in enumerate(page_rows):
        header = ''
        header += ' ' * label_width
        header += Monospace.get_ruler_line(i*width + 1, i*width + len(page_row[0]))
        yield header
        for label, element_list in zip(labels, page_row):
            justified_label = Monospace.left_justify(label, label_width, ' ')
            yield ''.join([justified_label] + list(element_list))
        if i < len(page_rows) - 1:
            yield ''

def do_analysis(mixture_model, alignment, tree):
    """
    @param mixture_model: a mixture of nucleotide rate matrices
    @param alignment: a nucleotide alignment
    @param tree: the phylogenetic tree with branch lengths
    @return: an html string representing a whole html file
    """
    # For each column of the alignment get the likelihood for each category.
    # The rest of the analysis can proceed from this data alone.
    likelihood_columns = []
    # create a hash table to help decorate the tree
    header_to_node = {}
    for header in alignment.headers:
        try:
            node = tree.get_unique_node(header)
        except Newick.NewickSearchError, e:
            raise HandlingError(e)
        header_to_node[header] = node
    # get the information for each column
    for column in alignment.columns:
        # decorate the tree with the ordered states of the current column
        for header, state in zip(alignment.headers, column):
            header_to_node[header].state = state
        # get the likelihood for each category
        likelihoods = []
        for p, matrix in zip(mixture_model.mixture_parameters, mixture_model.rate_matrices):
            likelihoods.append(p * matrix.get_likelihood(tree))
        likelihood_columns.append(likelihoods)
    # the likelihood_columns variable has everything we need to write the response
    # define the likelihood legend
    likelihood_column_sums = [sum(likelihoods) for likelihoods in likelihood_columns]
    likelihood_legend = HeatMap.Legend(likelihood_column_sums, 5, 'L', HeatMap.white_red_gradient)
    # get the mixture for each column implied by the likelihoods at the column
    mixture_columns = []
    for likelihoods in likelihood_columns:
        total = sum(likelihoods)
        mixtures = [likelihood / total for likelihood in likelihoods]
        mixture_columns.append(mixtures)
    # get the conditional mixtures for the whole alignment
    total_mixture = []
    for proportions in zip(*mixture_columns):
        total_mixture.append(sum(proportions) / len(alignment.columns))
    # define the mixture legend
    mixture_legend = HeatMap.Legend(Util.flattened_nonrecursive(mixture_columns), 5, 'M', HeatMap.white_blue_gradient)
    # start writing the web page
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<head>'
    print >> out, '<style>'
    for legend in (likelihood_legend, mixture_legend):
        for line in legend.gen_style_lines():
            print >> out, line
    print >> out, '</style>'
    print >> out, '</head>'
    print >> out, '<body>'
    # write the log likelihood
    log_likelihood = sum(math.log(sum(likelihoods)) for likelihoods in likelihood_columns)
    print >> out, 'log likelihood:'
    print >> out, '<br/>'
    print >> out, '%f' % log_likelihood
    # write the log likelihood per column
    print >> out, '<br/><br/>'
    print >> out, 'log likelihood per column:'
    print >> out, '<br/>'
    print >> out, '%f' % (log_likelihood / len(alignment.columns))
    # write the conditional mixtures for the whole alignment
    print >> out, '<br/><br/>'
    print >> out, 'conditional mixture:'
    print >> out, '<br/>'
    for proportion in total_mixture:
        print >> out, '%f</br>' % proportion
    # begin the pre environment
    print >> out, '<pre>'
    # write the alignment
    labels = alignment.headers + ['category 1', 'category 2', 'likelihood']
    element_lists = [list(seq) for seq in alignment.sequences]
    for proportions in zip(*mixture_columns):
        mixture_elements = []
        for proportion in proportions:
            mixture_elements.append('<span class="%s"> </span>' % mixture_legend.value_to_css_class(proportion))
        element_lists.append(mixture_elements)
    likelihood_elements = []
    for likelihood in likelihood_column_sums:
        likelihood_elements.append('<span class="%s"> </span>' % likelihood_legend.value_to_css_class(likelihood))
    element_lists.append(likelihood_elements)
    for line in do_analysis_helper(labels, element_lists, 60):
        print >> out, line
    # write the legend
    print >> out, ''
    print >> out, 'mixture key:'
    for line in reversed(list(mixture_legend.gen_legend_lines())):
        print >> out, line
    print >> out, ''
    print >> out, 'likelihood key:'
    for line in reversed(list(likelihood_legend.gen_legend_lines())):
        print >> out, line
    # end the pre environment
    print >> out, '</pre>'
    # terminate the file
    print >> out, '</body>'
    print >> out, '</html>'
    return out.getvalue()

def get_kappa(index):
    return [2, 2][index]

def get_weight(index):
    return [3, 1][index]

def get_frequency_string(index):
    category_to_frequencies = [
            [1, 1, 1, 1],
            [1, 4, 4, 1]
            ]
    frequencies = category_to_frequencies[index]
    arr = []
    for nt, frequency in zip('ACGT', frequencies):
        arr.append('%s : %s' % (nt, frequency))
    return '\n'.join(arr)
