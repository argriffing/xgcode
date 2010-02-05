"""Evaluate a codon alignment given a tree and a Direct Protein mixture model.

Find the probability of a codon alignment
given a tree and a Direct Protein mixture model.
The branch length is the expected number of codon substitutions on a branch.
See "Population Genetics Without Intraspecific Data" by Thorne et al.
for more information about the Direct Protein model.
"""

import math
from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Discretizer
import Util
import Monospace
import HeatMap
import Fasta
import DirectProtein
import Form
import iterutils

# This codon alignment is unrealistic.
g_sample_alignment_string = """
>Gibbon
TGTGAAGTGCTGCCTTTATTCTACTCCCCTTACCCTAGATTCATTCGTTTAGAGGGTGCC
CAGATCAGACGGTTCATAATAGGCCCGGAAGGCCACGGAGGACGATATCCACCCCTCTTC
CTATGCAATTCATGTCCTAGGGTCTCTTAT
>Orangutan
TGTGAAGTGGTGCCTCTATTCTACACCCCTTATCGTCGATTCATTAGTCTAGAAGGTGCC
CAAATTAGACGGTTCATAATAGGCCCGCGAGGGTGTGCAGGTCGATTTCCACCCCTATAT
TCGTGCAATTTCTCCCCTAGGGTCTCCTAT
>Gorilla
TATGATGTACCGCCTTTACTCTACATCCCACACGCGAAGCTCATTAGTCCAGAAGAGTCC
CAGAACAGAGTGGTCGTAAGGACCGCTGGAAGCCATGCAGGTAGATTCCCGCCCTTTACC
CCGTGCAATGTCCTCCCTATTGTTTCATAT
>Chimpanzee
TGCGAGAGGGCTCGTTTACTGTACATGCCTTGCTCTCAACCTACTAGTCCAGGGGGGCCC
TTTAACACATCTACCGTGATGAGCTCCGTATGGCATGCAGGTAGATTCCCGCTCTTAACC
TCTTGTAATGTCTTCCGTTTTGTTTCATAT
>Human
TGCGAGAGGGCTCCTTTACTGTACATGCCTTACTCTCAACCCATTAGTCCAGAGGGGCCC
CTTAGCAAACTTACCGTGATGAGCTCCGTTTGGCATGCAGGTAGATTCCCGCCCTTAACC
TCCTGTAATGTCTTCCGTTTCGTTTCATAT
"""

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
            Form.MultiLine('alignment', 'codon alignment',
                g_sample_alignment_string.strip()),
            Form.MultiLine('model', 'Direct Protein mixture model',
                DirectProtein.get_sample_xml_string().strip())]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the nucleotide alignment
    try:
        alignment = Fasta.CodonAlignment(StringIO(fs.alignment))
    except Fasta.AlignmentError, e:
        raise HandlingError(e)
    # get the normalized Direct Protein mixture model
    mixture_model = DirectProtein.deserialize_mixture_model(fs.model)
    mixture_model.normalize()
    # get the html string
    html_string = do_analysis(mixture_model, alignment, tree)
    # write the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, html_string

def do_analysis_helper(labels, element_lists, w):
    """
    Chop up the rows of data.
    Yield lines of text to be displayed in an html pre tag.
    @param labels: row labels to be left justified
    @param element_lists: each row where each element is a codon or a span
    @param w: the width; the number of elements allowed per page row
    """
    if len(set(len(element_list) for element_list in element_lists)) != 1:
        msg = 'each element list should have the same nonzero length'
        raise ValueError(msg)
    label_width = max(len(label) for label in labels) + 1
    chopped_element_lists = [list(iterutils.chopped(element_list, w))
            for element_list in element_lists]
    page_rows = zip(*chopped_element_lists)
    for i, page_row in enumerate(page_rows):
        header = ''
        header += ' ' * label_width
        k = len(page_row[0])
        header += Monospace.get_codon_ruler_line(i*w + 1, i*w + k)
        yield header
        for label, element_list in zip(labels, page_row):
            justified_label = Monospace.left_justify(label, label_width, ' ')
            yield ''.join([justified_label] + list(element_list))
        if i < len(page_rows) - 1:
            yield ''

def do_analysis(mixture_model, alignment, tree):
    """
    @param mixture_model: a mixture of codon rate matrices
    @param alignment: a codon nucleotide alignment
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
        for p, matrix in zip(mixture_model.mixture_parameters,
                mixture_model.rate_matrices):
            likelihoods.append(p * matrix.get_likelihood(tree))
        likelihood_columns.append(likelihoods)
    # The likelihood_columns variable
    # has everything we need to write the response.
    # Define the likelihood legend.
    likelihood_column_sums = [sum(likelihoods)
            for likelihoods in likelihood_columns]
    likelihood_legend = HeatMap.Legend(likelihood_column_sums,
            5, 'L', HeatMap.white_red_gradient)
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
    flattened_columns = Util.flattened_nonrecursive(mixture_columns)
    mixture_legend = HeatMap.Legend(flattened_columns,
            5, 'M', HeatMap.white_blue_gradient)
    # start writing the web page
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<head>'
    print >> out, '<style type="text/css">'
    print >> out, '.oddcodon {background-color: rgb(255, 255, 200)}'
    for legend in (likelihood_legend, mixture_legend):
        for line in legend.gen_style_lines():
            print >> out, line
    print >> out, '</style>'
    print >> out, '</head>'
    print >> out, '<body>'
    # write the log likelihood
    log_likelihood = sum(math.log(sum(likelihoods))
            for likelihoods in likelihood_columns)
    print >> out, 'log likelihood:'
    print >> out, '<br/>'
    print >> out, '%f' % log_likelihood
    # write the log likelihood per column
    print >> out, '<br/><br/>'
    print >> out, 'log likelihood per codon column:'
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
    labels = alignment.headers
    labels += ['category 1', 'category 2', 'category 3', 'likelihood']
    element_lists = []
    for codons in alignment.sequences:
        element_list = []
        for i, codon in enumerate(codons):
            if i % 2:
                element = '<span class="oddcodon">' + codon + '</span>'
            else:
                element = codon
            element_list.append(element)
        element_lists.append(element_list)
    for proportions in zip(*mixture_columns):
        mixture_elements = []
        for proportion in proportions:
            css_class = mixture_legend.value_to_css_class(proportion)
            mixture_elements.append('<span class="%s">   </span>' % css_class)
        element_lists.append(mixture_elements)
    likelihood_elements = []
    for likelihood in likelihood_column_sums:
        css_class = likelihood_legend.value_to_css_class(likelihood)
        likelihood_elements.append('<span class="%s">   </span>' % css_class)
    element_lists.append(likelihood_elements)
    for line in do_analysis_helper(labels, element_lists, 20):
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
