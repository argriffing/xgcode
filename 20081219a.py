"""Get intermediate values from a single column MAPP calculation.
"""

import cgi
from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import FelTree
import Newick
import NewickIO
import MatrixUtil
import HtmlTable
import MAPP
import Codon
import LeafWeights
import const

g_default_string = const.read('20100730s')

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree
    default_tree = NewickIO.parse(g_default_string, FelTree.NewickTree)
    default_tree_string = NewickIO.get_narrow_newick_string(default_tree, 60)
    # Define the default lines
    # that tell the amino acids in the column of the alignment.
    default_alignment_lines = [
            'Hg18 F',
            'bosTau3 V',
            'canFam2 V',
            'danRer5 L',
            'galGal3 A',
            'mm9 L',
            'monDom4 H',
            'panTro2 F',
            'rheMac2 F',
            'rn4 L']
    default_alignment_string = '\n'.join(default_alignment_lines)
    # define the list of form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', default_tree_string),
            Form.MultiLine('column',
                'amino acid of each taxon', default_alignment_string),
            Form.CheckGroup('options',
                'show these intermediate values', [
                Form.CheckItem('show_raw_pc_table',
                    'the raw physicochemical property table', True),
                Form.CheckItem('show_standardized_pc_table',
                    'the standardized physicochemical property table', True),
                Form.CheckItem('show_pc_correlation_matrix',
                    'the physicochemical property correlation matrix', True),
                Form.CheckItem('show_tree',
                    'the pruned phylogenetic tree', True),
                Form.CheckItem('show_weights',
                    'the taxon weights', True),
                Form.CheckItem('show_aa_distribution',
                    'the estimated amino acid distribution', True),
                Form.CheckItem('show_pc_distribution',
                    'the estimated physicochemical property distn', True),
                Form.CheckItem('show_deviations',
                    'the aa physicochemical property deviations', True),
                Form.CheckItem('show_impact_scores',
                    'the impact score for each amino acid', True),
                Form.CheckItem('show_p_values',
                    'the p-value for each amino acid', True)])]
    return form_objects

def get_form_out():
    return FormOut.Html()

def aa_letter_to_aa_index(aa_letter):
    """
    @param aa_letter: an amino acid letter
    @return: None or an index between 0 and 19
    """
    for i, aa in enumerate(Codon.g_aa_letters):
        if aa == aa_letter:
            return i
    return None

def get_tree_and_column(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: the pruned tree and a map from taxa to amino acids
    """
    # get the newick tree.
    tree = NewickIO.parse(fs.tree, Newick.NewickTree)
    unordered_tip_names = set(node.name for node in tree.gen_tips())
    # get the lines that give an amino acid for each of several taxa
    column_lines = Util.get_stripped_lines(StringIO(fs.column))
    if len(column_lines) < 7:
        msg = 'the alignment column should have at least seven taxa'
        raise HandlingError(msg)
    # get the mapping from taxon to amino acid
    taxon_to_aa_letter = {}
    for line in column_lines:
        pair = line.split()
        if len(pair) != 2:
            raise HandlingError('invalid line: %s' % line)
        taxon, aa_letter = pair
        aa_letter = aa_letter.upper()
        if aa_letter not in Codon.g_aa_letters:
            msg = 'expected an amino acid instead of this: %s' % aa_letter
            raise HandlingError(msg)
        taxon_to_aa_letter[taxon] = aa_letter
    # Assert that the names in the column are a subset of the names
    # of the tips of the tree.
    unordered_taxon_names = set(taxon_to_aa_letter)
    weird_names = unordered_taxon_names - unordered_tip_names
    if weird_names:
        msg = 'these taxa were not found on the tree: %s' % str(weird_names)
        raise HandlingError(msg)
    # define the taxa that will be pruned
    names_to_remove = unordered_tip_names - unordered_taxon_names
    # prune the tree
    for name in names_to_remove:
        tree.prune(tree.get_unique_node(name))
    # merge segmented branches 
    internal_nodes_to_remove = [node for node in tree.preorder()
            if node.get_child_count() == 1] 
    for node in internal_nodes_to_remove: 
        tree.remove_node(node) 
    return tree, taxon_to_aa_letter

def get_response_content(fs):
    # start writing the html response
    out = StringIO()
    print >> out, '<html>'
    print >> out, '<body>'
    # get the tree and the column sent by the user
    pruned_tree, taxon_to_aa_letter = get_tree_and_column(fs)
    # get the weights of the taxa
    taxon_weight_pairs = LeafWeights.get_stone_weights(pruned_tree)
    # show the raw physicochemical property table
    if fs.show_raw_pc_table:
        print >> out, 'raw physicochemical property table:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = Codon.g_aa_letters
        table = MAPP.g_property_array
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the standardized physicochemical property table
    standardized_property_array = MAPP.get_standardized_property_array(
            MAPP.g_property_array)
    # show the standardized physicochemical property table
    if fs.show_standardized_pc_table:
        print >> out, 'standardized physicochemical property table:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = Codon.g_aa_letters
        table = standardized_property_array
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the physicochemical property correlation matrix
    correlation_matrix = MAPP.get_property_correlation_matrix(
            standardized_property_array)
    # show the physicochemical property correlation matrix
    if fs.show_pc_correlation_matrix:
        print >> out, 'physicochemical property correlation matrix:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = MAPP.g_property_names
        table = correlation_matrix
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # show the pruned tree
    if fs.show_tree:
        tree_string = NewickIO.get_narrow_newick_string(pruned_tree, 80)
        lines = StringIO(tree_string).readlines()
        lines = [line.rstrip() for line in lines]
        print >> out, 'pruned phylogenetic tree in newick format:'
        print >> out, '<pre>'
        for line in lines:
            print >> out, cgi.escape(line)
        print >> out, '</pre>'
        print >> out, '<br/>'
    # show the weights
    if fs.show_weights:
        taxa, weights = zip(*taxon_weight_pairs)
        table = [weights]
        row_labels = ['weight']
        col_labels = taxa
        print >> out, 'taxon weights:'
        print >> out, '<br/>'
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # estimate the amino acid distribution for the column,
    # taking into account the tree and a uniform prior.
    weights = []
    aa_indices = []
    for taxon, weight in taxon_weight_pairs:
        weights.append(weight)
        aa_indices.append(aa_letter_to_aa_index(taxon_to_aa_letter[taxon]))
    aa_distribution = MAPP.estimate_aa_distribution(weights, aa_indices)
    # show the estimated amino acid distribution
    if fs.show_aa_distribution:
        table = [aa_distribution]
        row_labels = ['weight']
        col_labels = Codon.g_aa_letters
        print >> out, 'estimated amino acid distribution:'
        print >> out, '<br/>'
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # estimate the mean and variance of each physicochemical property
    est_pc_means = MAPP.estimate_property_means(
            standardized_property_array, aa_distribution)
    est_pc_variances = MAPP.estimate_property_variances(
            standardized_property_array, aa_distribution)
    # show the estimated mean and variance of each physicochemical property
    if fs.show_pc_distribution:
        table = [est_pc_means, est_pc_variances]
        row_labels = ['mean', 'variance']
        col_labels = MAPP.g_property_names
        print >> out, 'estimated physicochemical property moments:'
        print >> out, '<br/>'
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the deviation from each property mean
    # for each possible amino acid
    deviations = MAPP.get_deviations(
            est_pc_means, est_pc_variances, standardized_property_array)
    # show the deviation from each property mean for each possible amino acid
    if fs.show_deviations:
        print >> out, 'deviations of amino acids from the normal distribution'
        print >> out, 'estimated for each property:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = Codon.g_aa_letters
        table = deviations
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the impact scores
    impact_scores = MAPP.get_impact_scores(correlation_matrix, deviations)
    # show the impact scores
    if fs.show_impact_scores:
        table = [impact_scores]
        row_labels = ['impact']
        col_labels = Codon.g_aa_letters
        print >> out, 'impact scores:'
        print >> out, '<br/>'
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the p-values
    p_values = []
    for score in impact_scores:
        ntaxa = len(taxon_weight_pairs)
        p_values.append(MAPP.get_p_value(score, ntaxa))
    # show the p-values
    if fs.show_p_values:
        table = [p_values]
        row_labels = ['p-value']
        col_labels = Codon.g_aa_letters
        print >> out, 'p-values:'
        print >> out, '<br/>'
        print >> out, HtmlTable.get_labeled_table_string(
                col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # write the html footer
    print >> out, '</body>'
    print >> out, '</html>'
    # return the response
    return out.getvalue()
