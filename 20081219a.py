"""Get intermediate values from a single column MAPP calculation.
"""

import cgi
from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import Util
import FelTree
import Newick
import NewickIO
import MatrixUtil
import HtmlTable
import MAPP
import Codon
import LeafWeights

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree
    default_tree_lines = [
            '((((((((',
            '(((mm9:0.076274,rn4:0.084383):0.200607,',
            'cavPor2:0.202990):0.034350,',
            'oryCun1:0.208548):0.014587,',
            '((((((Hg18:0.005873,panTro2:0.007668):0.013037,',
            'ponAbe2:0.02):0.013037,rheMac2:0.031973):0.0365,',
            'calJac1:0.07):0.0365,otoGar1:0.151185):0.015682,',
            'tupBel1:0.162844):0.006272):0.019763,',
            '((sorAra1:0.248532,eriEur1:0.222255):0.045693,',
            '(((canFam2:0.101137,felCat3:0.098203):0.048213,',
            'equCab1:0.099323):0.007287,',
            'bosTau3:0.163945):0.012398):0.018928):0.030081,',
            '(dasNov1:0.133274,(loxAfr1:0.103030,',
            'echTel1:0.232706):0.049511):0.008424):0.213469,',
            'monDom4:0.320721):0.088647,',
            'ornAna1:0.488110):0.118797,',
            '(galGal3:0.395136,anoCar1:0.513962):0.093688):0.151358,',
            'xenTro2:0.778272):0.174596,',
            '(((tetNig1:0.203933,fr1:0.239587):0.203949,',
            '(gasAcu1:0.314162,oryLat1:0.501915):0.055354):0.346008,',
            'danRer5:0.730028):0.174596);']
    tree_string = '\n'.join(default_tree_lines)
    default_tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    default_tree_string = NewickIO.get_narrow_newick_string(default_tree, 60)
    # define the default lines that tell the amino acids in the column of the alignment
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
            Form.MultiLine('column', 'amino acid of each taxon', default_alignment_string),
            Form.CheckGroup('options', 'show these intermediate values', [
                Form.CheckItem('show_raw_pc_table', 'the raw physicochemical property table', True),
                Form.CheckItem('show_standardized_pc_table', 'the standardized physicochemical property table', True),
                Form.CheckItem('show_pc_correlation_matrix', 'the physicochemical property correlation matrix', True),
                Form.CheckItem('show_tree', 'the pruned phylogenetic tree', True),
                Form.CheckItem('show_weights', 'the taxon weights', True),
                Form.CheckItem('show_aa_distribution', 'the estimated amino acid distribution', True),
                Form.CheckItem('show_pc_distribution', 'the estimated physicochemical property distribution', True),
                Form.CheckItem('show_deviations', 'the amino acid physicochemical property deviations', True),
                Form.CheckItem('show_impact_scores', 'the impact score for each amino acid', True),
                Form.CheckItem('show_p_values', 'the p-value for each amino acid', True)])]
    return form_objects

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

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
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
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the standardized physicochemical property table
    standardized_property_array = MAPP.get_standardized_property_array(MAPP.g_property_array)
    # show the standardized physicochemical property table
    if fs.show_standardized_pc_table:
        print >> out, 'standardized physicochemical property table:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = Codon.g_aa_letters
        table = standardized_property_array
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the physicochemical property correlation matrix
    correlation_matrix = MAPP.get_property_correlation_matrix(standardized_property_array)
    # show the physicochemical property correlation matrix
    if fs.show_pc_correlation_matrix:
        print >> out, 'physicochemical property correlation matrix:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = MAPP.g_property_names
        table = correlation_matrix
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
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
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
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
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # estimate the mean and variance of each physicochemical property
    est_pc_means = MAPP.estimate_property_means(standardized_property_array, aa_distribution)
    est_pc_variances = MAPP.estimate_property_variances(standardized_property_array, aa_distribution)
    # show the estimated mean and variance of each physicochemical property
    if fs.show_pc_distribution:
        table = [est_pc_means, est_pc_variances]
        row_labels = ['mean', 'variance']
        col_labels = MAPP.g_property_names
        print >> out, 'estimated physicochemical property moments:'
        print >> out, '<br/>'
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # calculate the deviation from each property mean for each possible amino acid
    deviations = MAPP.get_deviations(est_pc_means, est_pc_variances, standardized_property_array)
    # show the deviation from each property mean for each possible amino acid
    if fs.show_deviations:
        print >> out, 'deviations of amino acids from the normal distribution estimated for each property:'
        print >> out, '<br/>'
        col_labels = MAPP.g_property_names
        row_labels = Codon.g_aa_letters
        table = deviations
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
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
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
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
        print >> out, HtmlTable.get_labeled_table_string(col_labels, row_labels, table)
        print >> out, '<br/><br/>'
    # write the html footer
    print >> out, '</body>'
    print >> out, '</html>'
    # write the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, out.getvalue()

