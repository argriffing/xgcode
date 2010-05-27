"""Do a MAPP analysis.
"""

from StringIO import StringIO
import itertools
import random
import time
import sys

import argparse
import numpy as np

from SnippetUtil import HandlingError
import Form
import Util

g_default_string = """
foo 1 1 1
bar 1 1 1
baz 1 0 1
""".strip()


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('words',
                'a list of OTUs names and binary character vectors',
                g_default_string),
            Form.Integer('nwords',
                'find a subset with this many OTUs', 2),
            Form.Integer('nchars',
                'find a subset with this many binary characters', 1)]
    return form_objects


def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    lines = Util.get_stripped_lines(StringIO(fs.words))
    try:
        words = [Word(line) for line in lines]
        validate_words(words)
    except WordError, e:
        raise HandlingError(e)
    text = process(words, fs.nwords, fs.nchars)
    return [('Content-Type', 'text/plain')], text


class WordError(Exception): pass

class Word(object):
    def __init__(self, line):
        # store the original line
        self.line = line
        # extract the name and the binary vector
        elements = line.split()
        self.name = elements[0]
        for x in elements[1:]:
            if x not in ('0', '1'):
                msg = 'expected 0 or 1 but got ' + x
                raise WordError(msg)
        self.v = np.array([int(x) for x in elements[1:]])

def validate_words(words):
    """
    Make sure the binary vectors are the same lengths.
    @param words: word objects
    """
    lengths = set(len(w.v) for w in words)
    if len(lengths) != 1:
        msg = 'each binary vector should be the same length'
        raise WordError(msg)

def words_to_matrix(words):
    """
    @param words: validated word objects
    @return: matrix where rows are OTUs and cols are binary characters
    """
    nrows = len(words)
    ncols = len(words[0].v)
    M = np.zeros((nrows, ncols), int)
    for i, w in enumerate(words):
        for j, x in enumerate(w.v):
            M[i,j] = x
    return M

def get_separation(M, row_indices, col_indices):
    """
    @param M: a binary matrix
    @param row_indices: a set of selected row indices
    @param col_indices: a set of selected column indices
    @return: the min hamming distance among pairs of selected rows
    """
    # first extract some vectors from the matrix
    vectors = []
    for i in row_indices:
        v = np.array([M[i,j] for j in col_indices])
        vectors.append(v)
    # now get the distance
    pairs = itertools.combinations(vectors, 2)
    return min(np.dot(vb - va, vb - va) for va, vb in pairs)

def get_selections(M, nrows, ncols, nseconds):
    """
    Select a set of rows and a set of columns.
    The selections should give a large separation,
    which I am defining as the minimum hamming distance
    among all pairs of selected rows.
    @param M: a binary matrix; rows are OTUs and cols are binary characters
    @param nrows: the number of requested rows
    @param ncols: the number of requested columns
    @param nseconds: a time limit for the search
    @return: a row index set and column index set
    """
    nrows_total, ncols_total = M.shape
    best_selections = None
    best_separation = None
    t = time.time()
    while time.time() - t < nseconds:
        row_indices = set(random.sample(range(nrows_total), nrows))
        col_indices = set(random.sample(range(ncols_total), ncols))
        d = get_separation(M, row_indices, col_indices)
        if (best_separation is None) or (best_separation < d):
            best_separation = d
            best_selections = (row_indices, col_indices)
    return best_selections

def process(words, nwords, nchars, nseconds=2.0):
    """
    @param words: a sequence of word objects
    @param nwords: find a subset of this many word objects
    @param nchars: find a subset of this many binary characters
    @param nseconds: a time limit
    @return: multiline string
    """
    out = StringIO()
    if len(words) < nwords:
        msg = 'the number of OTUs is smaller than the desired sample'
        raise HandlingError(msg)
    if len(words[0].v) < nchars:
        msg = 'the number of characters is smaller than the desired sample'
        raise HandlingError(msg)
    # create the matrix
    M = words_to_matrix(words)
    # select row and column indices
    row_indices, col_indices = get_selections(M, nwords, nchars, nseconds)
    sorted_row_indices = list(sorted(row_indices))
    sorted_col_indices = list(sorted(col_indices))
    # print the separation
    d = get_separation(M, row_indices, col_indices)
    print >> out, 'best separation:', d
    # print the index selections
    print >> out, 'selected row indices:', sorted_row_indices
    print >> out, 'selected column indices:', sorted_col_indices
    # print some selected values
    for i in sorted_row_indices:
        w = words[i]
        s = ' '.join(str(M[i,j]) for j in sorted_col_indices)
        print >> out, w.name + '\t' + s
    return out.getvalue().rstrip()

def process(raw_lines):
    out = StringIO()
    line = Util.get_first(Util.stripped_lines(raw_lines))
    otu_name, genotype_string = line.split(None, 1)
    genotypes = genotype_string.split()
    for i, genotype in enumerate(genotypes):
        name = 'SNP_' + str(i)
        chromosome = '1'
        morgans = '0.0'
        bases = i+1
        row = [name, chromosome, morgans, bases]
        print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue().rstrip()

def main(args):
    print process(sys.stdin)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    args = parser.parse_args()
    main(args)



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


class Snp(object):

    def __init__(self, snp_lines):
        # get the codon
        pass

    def parse_first_line(self, line):
        """
        @param line: a line of comma separated values
        """
        elements = line.split(',')
        elements = [x.strip(string.whitespace + '"') for x in elements]


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
    # calculate the standardized physicochemical property table
    standardized_property_array = MAPP.get_standardized_property_array(
            MAPP.g_property_array)
    # calculate the physicochemical property correlation matrix
    correlation_matrix = MAPP.get_property_correlation_matrix(
            standardized_property_array)
    # estimate the amino acid distribution for the column,
    # taking into account the tree and a uniform prior.
    weights = []
    aa_indices = []
    for taxon, weight in taxon_weight_pairs:
        weights.append(weight)
        aa_indices.append(aa_letter_to_aa_index(taxon_to_aa_letter[taxon]))
    aa_distribution = MAPP.estimate_aa_distribution(weights, aa_indices)
    # estimate the mean and variance of each physicochemical property
    est_pc_means = MAPP.estimate_property_means(
            standardized_property_array, aa_distribution)
    est_pc_variances = MAPP.estimate_property_variances(
            standardized_property_array, aa_distribution)
    # calculate the deviation from each property mean
    # for each possible amino acid
    deviations = MAPP.get_deviations(
            est_pc_means, est_pc_variances, standardized_property_array)
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
    # write the response
    response_headers = [('Content-Type', 'text/html')]
    return response_headers, out.getvalue()


