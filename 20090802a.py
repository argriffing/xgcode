"""Examine deep splits of the tree of life.

This is mostly for a figure for a paper.
The default data is from a project called the interactive tree of life.
It currently supports the Felsenstein-pruning-like way of getting
subsequent matrices given a split,
and now I'm working on support for a neighbor-joining-like
way of getting the subsequent distance matrices.
This latter approach is probably worse,
but it is easer to present because it
can be explained in parallel with neighbor joining.
"""

from StringIO import StringIO
import os
import zipfile

import numpy as np

from SnippetUtil import HandlingError
import NewickIO
import FelTree
import Newick
import BuildTreeTopology
import SchurAlgebra
import Euclid
import MatrixUtil
from Form import CheckGroup
from Form import CheckItem
from Form import RadioGroup
from Form import RadioItem
import Form
import FormOut
import const

g_archaea_data = const.read('20090802a')
g_bacteria_data = const.read('20090802b')
g_eukaryota_data = const.read('20090802c')
g_full_data = const.read('20090802d')

g_valid_domains = set(['archaea', 'bacteria', 'eukaryota'])

#FIXME multiple output types

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            RadioGroup('distance_options', 'recursive matrix construction', [
                RadioItem('like_pruning',
                    'go through the Laplacian, like Felsenstein pruning',
                    True),
                RadioItem('like_nj',
                    'directly use distances, like neighbor joining')]),
            CheckGroup('format_options', 'options', [
                CheckItem('supplementary',
                    'download the supplementary data', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def remove_redundant_nodes(tree):
    """
    @param tree: a Newick NewickTree
    """
    # find nodes with a branch length of zero
    zero_blen_nodes = [node for node in tree.preorder() if node.blen == 0.0]
    # deal with these nodes appropriately
    for node in zero_blen_nodes:
        if node.blen != 0.0:
            # it is possible that this node no longer has branch length zero
            continue
        elif node.has_children():
            # remove internal nodes with zero branch length
            tree.remove_node(node)
        else:
            # prune terminal nodes with zero branch length
            intersection = tree.prune(node)
            # remove intersection nodes that have a single child
            if intersection.get_child_count() == 1:
                tree.remove_node(intersection)

class SupplementaryObject:
    """
    An object of this class defines supplementary data for a paper.
    It will comprise three files that are zipped together.
    """

    def __init__(self, taxon_to_domain, newick_string, use_generalized_nj):
        """
        @param taxon_to_domain: maps a taxon name to a taxonomic classification
        @param newick_string: a string that defines a newick tree
        @param use_generalized_nj: True if we use an old method of outgrouping
        """
        # store the info that defines the supplementary data
        self.taxon_to_domain = taxon_to_domain
        self.newick_string = newick_string
        # create the full and pruned trees
        self.full_tree = None
        self.pruned_tree = None
        self._create_trees()
        # create the ordered names of taxa in the pruned tree
        self.pruned_names = None
        self._create_ordered_names()
        # create some supplementary objects
        self.first_split_object = None
        self.second_split_object = None
        self._do_analysis(use_generalized_nj)

    def _create_trees(self):
        """
        Create the full tree and the pruned tree.
        The full tree is a Newick.NewickTree,
        and the pruned tree is a FelTree.NewickTree object.
        """
        # create the full tree
        self.full_tree = NewickIO.parse(self.newick_string, Newick.NewickTree)
        # create the pruned tree through a temporary tree that will be modified
        temp_tree = NewickIO.parse(self.newick_string, Newick.NewickTree)
        remove_redundant_nodes(temp_tree)
        pruned_newick_string = NewickIO.get_newick_string(temp_tree)
        self.pruned_tree = NewickIO.parse(pruned_newick_string, FelTree.NewickTree)

    def _create_ordered_names(self):
        """
        Create the ordered list of names.
        """
        pruned_archaea_names = self._get_pruned_names_in_domain('archaea')
        pruned_bacteria_names = self._get_pruned_names_in_domain('bacteria')
        pruned_eukaryota_names = self._get_pruned_names_in_domain('eukaryota')
        self.pruned_names = pruned_bacteria_names + pruned_archaea_names + pruned_eukaryota_names

    def _get_pruned_names_in_domain(self, domain):
        """
        @param domain: get names from this domain from the pruned tree
        @return: a list of names
        """
        assert domain in g_valid_domains
        unordered_pruned_names = list(node.get_name() for node in self.pruned_tree.gen_tips())
        return [x for x in unordered_pruned_names if self.taxon_to_domain[x] == domain]

    def _get_full_names_in_domain(self, domain):
        """
        @param domain: get names from this domain from the full tree
        @return: a list of names
        """
        assert domain in g_valid_domains
        unordered_names = list(node.get_name() for node in self.full_tree.gen_tips())
        return [x for x in unordered_names if self.taxon_to_domain[x] == domain]

    def _get_domains(self, names):
        """
        @param names: a collection of taxon names
        @return: the set of domains represented by these names
        """
        return set(self.taxon_to_domain[name] for name in names)

    def _do_analysis(self, use_generalized_nj):
        """
        Do some splits of the tree.
        @param use_generalized_nj: True if we use an old method of outgrouping
        """
        # define the distance matrix
        D = np.array(self.pruned_tree.get_distance_matrix(self.pruned_names))
        # get the primary split of the criterion matrix
        L = Euclid.edm_to_laplacian(D)
        v = BuildTreeTopology.laplacian_to_fiedler(L)
        eigensplit = BuildTreeTopology.eigenvector_to_split(v)
        # assert that the first split cleanly separates the bacteria from the rest
        left_indices, right_indices = eigensplit
        left_domains = self._get_domains([self.pruned_names[x] for x in left_indices])
        right_domains = self._get_domains([self.pruned_names[x] for x in right_indices])
        if ('bacteria' in left_domains) and ('bacteria' in right_domains):
            raise HandlingError('bacteria were not defined by the first split')
        # now we have enough info to define the first supplementary csv file
        self.first_split_object = SupplementarySpreadsheetObject(self.pruned_names, L, v)
        # define the bacteria indices vs the non-bacteria indices for the second split
        if 'bacteria' in left_domains:
            bacteria_indices = left_indices
            non_bacteria_indices = right_indices
        elif 'bacteria' in right_domains:
            bacteria_indices = right_indices
            non_bacteria_indices = left_indices
        # get the secondary split of interest
        if use_generalized_nj:
            D_secondary = BuildTreeTopology.update_generalized_nj(D, bacteria_indices)
            L_secondary = Euclid.edm_to_laplacian(D_secondary)
        else:
            L_secondary = SchurAlgebra.mmerge(L, bacteria_indices)
        full_label_sets = [set([i]) for i in range(len(self.pruned_names))]
        next_label_sets = SchurAlgebra.vmerge(full_label_sets, bacteria_indices)
        v_secondary = BuildTreeTopology.laplacian_to_fiedler(L_secondary)
        eigensplit_secondary = BuildTreeTopology.eigenvector_to_split(v_secondary)
        left_subindices, right_subindices = eigensplit_secondary
        pruned_names_secondary = []
        for label_set in next_label_sets:
            if len(label_set) == 1:
                label = list(label_set)[0]
                pruned_names_secondary.append(self.pruned_names[label])
            else:
                pruned_names_secondary.append('all-bacteria')
        # assert that the second split cleanly separates the eukaryota from the rest
        left_subdomains = self._get_domains([pruned_names_secondary[x] for x in left_subindices])
        right_subdomains = self._get_domains([pruned_names_secondary[x] for x in right_subindices])
        if ('eukaryota' in left_subdomains) and ('eukaryota' in right_subdomains):
            raise HandlingError('eukaryota were not defined by the second split')
        # now we have enough info to define the second supplementary csv file
        self.second_split_object = SupplementarySpreadsheetObject(pruned_names_secondary, L_secondary, v_secondary)

    def _get_name_summary(self):
        """
        @return: a multiline string
        """
        out = StringIO()
        # create an interestingly ordered list of tips of the full tree
        full_archaea_names = self._get_full_names_in_domain('archaea')
        full_bacteria_names = self._get_full_names_in_domain('bacteria')
        full_eukaryota_names = self._get_full_names_in_domain('eukaryota')
        full_names = full_bacteria_names + full_archaea_names + full_eukaryota_names
        # create an interestingly ordered list of tips of the pruned tree
        pruned_archaea_names = self._get_pruned_names_in_domain('archaea')
        pruned_bacteria_names = self._get_pruned_names_in_domain('bacteria')
        pruned_eukaryota_names = self._get_pruned_names_in_domain('eukaryota')
        pruned_names = pruned_bacteria_names + pruned_archaea_names + pruned_eukaryota_names
        # show a summary of the original data
        print >> out, 'data summary before removing branches with zero length:'
        print >> out, len(full_archaea_names), 'archaea names in the original tree'
        print >> out, len(full_bacteria_names), 'bacteria names in the original tree'
        print >> out, len(full_eukaryota_names), 'eukaryota names in the original tree'
        print >> out, len(full_names), 'total names in the original tree'
        print >> out
        # show a summary of the processed data
        print >> out, 'data summary after removing branches with zero length:'
        print >> out, len(pruned_archaea_names), 'archaea names in the original tree'
        print >> out, len(pruned_bacteria_names), 'bacteria names in the original tree'
        print >> out, len(pruned_eukaryota_names), 'eukaryota names in the original tree'
        print >> out, len(pruned_names), 'total names in the processed non-degenerate tree'
        # return the multiline string
        return out.getvalue().strip()

    def get_verbose_summary(self):
        """
        @return: a multiline string
        """
        # begin the response
        out = StringIO()
        # show the number of taxa in various domains
        print >> out, self._get_name_summary()
        print >> out
        # show the pruned full tree
        formatted_tree_string = NewickIO.get_narrow_newick_string(self.pruned_tree, 120) 
        print >> out, 'this is the tree that represents all clades but for which redundant nodes have been pruned:'
        print >> out, formatted_tree_string
        print >> out
        # split the distance matrix
        D = np.array(self.pruned_tree.get_distance_matrix(self.pruned_names))
        L = Euclid.edm_to_laplacian(D)
        v = BuildTreeTopology.laplacian_to_fiedler(L)
        eigensplit = BuildTreeTopology.eigenvector_to_split(v)
        # report the eigendecomposition
        print >> out, get_eigendecomposition_report(D)
        print >> out
        # report the clade intersections of sides of the split
        side_names = [set(self.pruned_names[i] for i in side) for side in eigensplit]
        print >> out, 'domains represented by each side of the primary split:'
        print >> out, 'the left side has:\t', ', '.join(self._get_domains(side_names[0]))
        print >> out, 'the right side has:\t', ', '.join(self._get_domains(side_names[1]))
        print >> out
        # prepare to do the secondary splits
        left_indices, right_indices = eigensplit
        full_label_sets = [set([i]) for i in range(len(self.pruned_names))]
        # do the secondary splits
        for index_selection, index_complement in ((left_indices, right_indices), (right_indices, left_indices)):
            L_secondary = SchurAlgebra.mmerge(L, index_complement)
            next_label_sets = SchurAlgebra.vmerge(full_label_sets, index_complement)
            v = BuildTreeTopology.laplacian_to_fiedler(L_secondary)
            left_subindices, right_subindices = BuildTreeTopology.eigenvector_to_split(v)
            left_sublabels = set()
            for i in left_subindices:
                left_sublabels.update(next_label_sets[i])
            right_sublabels = set()
            for i in right_subindices:
                right_sublabels.update(next_label_sets[i])
            left_subnames = set(self.pruned_names[i] for i in left_sublabels)
            right_subnames = set(self.pruned_names[i] for i in right_sublabels)
            print >> out, 'domains represented by a subsplit:'
            print >> out, 'the left side has:\t', ', '.join(self._get_domains(left_subnames))
            print >> out, 'the right side has:\t', ', '.join(self._get_domains(right_subnames))
            print >> out
        # return the multiline string
        return out.getvalue().strip()

    def get_newick_file_contents(self):
        """
        @return: the contents of the supplementary newick file
        """
        return self.newick_string

    def get_first_split_file_contents(self):
        """
        @return: the contents of a supplementary csv file
        """
        return self.first_split_object.get_csv_contents(self.taxon_to_domain)

    def get_second_split_file_contents(self):
        """
        @return: the contents of a supplementary csv file
        """
        return self.second_split_object.get_csv_contents(self.taxon_to_domain)


class SupplementarySpreadsheetObject:
    """
    Each object of this class is associated with a spreadsheet.
    The supplementary data for my paper will include some spreadsheets.
    """

    def __init__(self, taxon_list, criterion_matrix, eigenvector):
        """
        @param taxon_list: an ordered list of taxon names
        @param criterion_matrix: a laplacian matrix as a numpy array
        @param eigenvector: the fiedler-like eigenvector for the criterion matrix
        """
        # do validation
        nrows, ncols = criterion_matrix.shape
        assert len(taxon_list) == nrows == ncols == len(eigenvector)
        # save the defining data
        self.taxon_list = taxon_list
        self.criterion_matrix = criterion_matrix
        self.eigenvector = list(eigenvector)

    def get_csv_contents(self, taxon_to_domain):
        """
        This member function creates the contents of a csv file.
        The first column is the list of taxon names.
        The second column is a phylogenetic classification of the taxon.
        The third column is an eigenvector of the criterion matrix.
        The remaining columns are columns of the criterion matrix.
        @param taxon_to_domain: maps a taxon name to a taxonomic classification
        @return: the contents of a csv file in convenient string form
        """
        # do validation
        assert set(taxon_to_domain.values()) <= set(['eukaryota', 'bacteria', 'archaea'])
        # define the columns of the csv file
        columns = []
        columns.append(self.taxon_list[:])
        columns.append([taxon_to_domain[taxon] for taxon in self.taxon_list])
        columns.append(self.eigenvector)
        columns.extend(self.criterion_matrix.T.tolist())
        # define the rows of the csv file
        rows = zip(*columns)
        # create the csv contents
        out = StringIO()
        for row in rows:
            print >> out, ', '.join(str(x) for x in row)
        csv_contents = out.getvalue()
        # return the csv contents
        return csv_contents


def get_eigendecomposition_report(D):
    """
    @param D: a distance matrix
    @return: a multi-line string
    """
    out = StringIO()
    # get some intermediate matrices and vectors
    L = Euclid.edm_to_laplacian(D)
    laplacian_fiedler = BuildTreeTopology.laplacian_to_fiedler(L)
    distance_fiedler = BuildTreeTopology.edm_to_fiedler(D)
    eigensplit = BuildTreeTopology.eigenvector_to_split(laplacian_fiedler)
    # report the two eigenvalue lists that should be the same
    HDH = MatrixUtil.double_centered(D)
    HSH = -0.5 * HDH
    w_distance, vt_distance = np.linalg.eigh(HSH)
    print >> out, 'the laplacian-derived and distance-derived eigenvalues:'
    w_laplacian, vt_laplacian = np.linalg.eigh(L)
    for a, b in zip(sorted(w_laplacian), sorted(w_distance)):
        print >> out, a, '\t', b
    print >> out
    # report the two fiedler vectors that should be the same
    print >> out, 'the laplacian-derived and distance-derived fiedler vectors:'
    for a, b in zip(laplacian_fiedler, distance_fiedler):
        print >> out, a, '\t', b
    return out.getvalue().strip()

def get_standard_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # begin the response
    out = StringIO()
    # show a summary of the original data
    print >> out, 'data summary before removing branches with zero length:'
    print >> out, len(archaea_names), 'archaea names in the original tree'
    print >> out, len(bacteria_names), 'bacteria names in the original tree'
    print >> out, len(eukaryota_names), 'eukaryota names in the original tree'
    print >> out, len(all_names), 'total names in the original tree'
    print >> out
    # get the pruned full tree
    pruned_full_tree = get_pruned_tree(full_tree)
    ordered_names = list(node.get_name() for node in pruned_full_tree.gen_tips())
    # show a summary of the processed data
    print >> out, 'data summary after removing branches with zero length:'
    print >> out, len(ordered_names), 'total names in the processed non-degenerate tree'
    print >> out
    # draw the pruned full tree
    print >> out, 'this is the tree that represents all clades but for which redundant nodes have been pruned:'
    formatted_tree_string = NewickIO.get_narrow_newick_string(pruned_full_tree, 120) 
    print >> out, formatted_tree_string
    print >> out
    # split the distance matrix
    D = np.array(pruned_full_tree.get_distance_matrix(ordered_names))
    L = Euclid.edm_to_laplacian(D)
    v = BuildTreeTopology.laplacian_to_fiedler(L)
    eigensplit = BuildTreeTopology.eigenvector_to_split(v)
    # report the eigendecomposition
    print >> out, get_eigendecomposition_report(D)
    # report the clade intersections of sides of the split
    side_names = [set(ordered_names[i] for i in side) for side in eigensplit]
    clade_name_pairs = ((archaea_names, 'archaea'), (bacteria_names, 'bacteria'), (eukaryota_names, 'eukaryota'))
    print >> out, 'clade intersections with each side of the split:'
    for side, side_name in zip(side_names, ('left', 'right')):
        for clade, clade_name in clade_name_pairs:
            if clade & side:
                print >> out, 'the', side_name, 'side intersects', clade_name
    print >> out
    # prepare to do the secondary splits
    left_indices, right_indices = eigensplit
    full_label_sets = [set([i]) for i in range(len(ordered_names))]
    # get a secondary split
    for index_selection, index_complement in ((left_indices, right_indices), (right_indices, left_indices)):
        L_s1 = SchurAlgebra.mmerge(L, index_complement)
        next_label_sets = SchurAlgebra.vmerge(full_label_sets, index_complement)
        v = BuildTreeTopology.laplacian_to_fiedler(L_s1)
        left_subindices, right_subindices = BuildTreeTopology.eigenvector_to_split(v)
        left_sublabels = set()
        for i in left_subindices:
            left_sublabels.update(next_label_sets[i])
        right_sublabels = set()
        for i in right_subindices:
            right_sublabels.update(next_label_sets[i])
        left_subnames = set(ordered_names[i] for i in left_sublabels)
        right_subnames = set(ordered_names[i] for i in right_sublabels)
        print >> out, 'clade intersections with a subsplit:'
        for clade, clade_name in clade_name_pairs:
            if clade & left_subnames:
                print >> out, 'the left side intersects', clade_name
        for clade, clade_name in clade_name_pairs:
            if clade & right_subnames:
                print >> out, 'the right side intersects', clade_name
        print >> out
    # show debug info
    print >> out, 'archaea names:'
    print >> out, '\n'.join(x for x in sorted(archaea_names))
    print >> out
    print >> out, 'bacteria names:'
    print >> out, '\n'.join(x for x in sorted(bacteria_names))
    print >> out
    print >> out, 'eukaryota names:'
    print >> out, '\n'.join(x for x in sorted(eukaryota_names))
    print >> out
    # return the response
    response_text = out.getvalue().strip()
    return [('Content-Type', 'text/plain')], response_text

def newick_string_to_tip_names(newick_string):
    """
    Get a set of names from a newick string.
    @param newick_string: the newick string that defines a tree
    @return: the set of names at the tips of the tree
    """
    tree = NewickIO.parse(newick_string, Newick.NewickTree)
    names = set(node.get_name() for node in tree.gen_tips())
    return names

def get_supplementary_object(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a object with all of the information for the supplementary data
    """
    # extract the name sets from the newick tree strings
    archaea_names = newick_string_to_tip_names(g_archaea_data)
    bacteria_names = newick_string_to_tip_names(g_bacteria_data)
    eukaryota_names = newick_string_to_tip_names(g_eukaryota_data)
    all_names = newick_string_to_tip_names(g_full_data)
    # validate the sets of names
    nfull = len(all_names)
    ndisjoint = len(archaea_names) + len(bacteria_names) + len(eukaryota_names)
    if ndisjoint != nfull:
        msg_a = 'there are %d taxa in the full tree ' % nfull
        msg_b = 'but %d taxa in its subtrees' % ndisjoint
        raise HandlingError(msg_a + msg_b)
    disjoint_union = archaea_names | bacteria_names | eukaryota_names
    if disjoint_union != all_names:
        msg_a = 'the set of taxa in the full tree '
        msg_b = 'is not the union of taxa in its subtrees'
        raise HandlingError(msg_a + msg_b)
    # create the map from taxon name to taxonomic category
    taxon_to_domain = {}
    for name in archaea_names:
        taxon_to_domain[name] = 'archaea'
    for name in bacteria_names:
        taxon_to_domain[name] = 'bacteria'
    for name in eukaryota_names:
        taxon_to_domain[name] = 'eukaryota'
    taxon_to_domain['all-bacteria'] = 'bacteria'
    # create the supplementary object
    use_generalized_nj = fs.like_nj
    supplementary_object = SupplementaryObject(
            taxon_to_domain, g_full_data, use_generalized_nj)
    # return the supplementary object
    return supplementary_object

def get_supplementary_response(supplementary_object):
    """
    @param supplementary_object: the three files can be readily obtained from this object
    @return: a (response_headers, response_text) pair
    """
    # create the zipfile
    fout = StringIO()
    zout = zipfile.ZipFile(fout, mode='w', compression=zipfile.ZIP_DEFLATED)
    zout.writestr('all.tree', supplementary_object.get_newick_file_contents())
    zout.writestr('first-split.csv', supplementary_object.get_first_split_file_contents())
    zout.writestr('second-split.csv', supplementary_object.get_second_split_file_contents())
    zout.close()
    response_text = fout.getvalue()
    # return the response
    response_headers = [
            ('Content-Type', 'application/zip'),
            ('Content-Disposition', 'attachment; filename=sup.zip')]
    return response_headers, response_text

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # do the analysis
    supplementary_object = get_supplementary_object(fs)
    # return the requested view of the first two interesting splits
    if fs.supplementary:
        return get_supplementary_response(supplementary_object)
    else:
        response_text = supplementary_object.get_verbose_summary()
        return [('Content-Type', 'text/plain')], response_text
