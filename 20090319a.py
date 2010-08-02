"""Calculate modulated network modularity given MMC input and output.

MMC is Modulated Modularity Clustering as in the paper
"Modulated Modularity Clustering as an Exploratory Tool
for Functional Genomic Inference".
"""

import math

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import Util
import MatrixUtil
import const

g_karate_data = const.read('20090406a')


def get_form():
    """
    @return: the body of a form
    """
    # define the default observation lines
    default_observation_matrix = [
            ['Gene/Observation', 'Observation 1', 'Observation 2', 'Observation 3', 'Observation 4', 'Observation 5'],
            ['Gene 1', '4.1', '6.23', '3.141529', '8.008', '1.2'],
            ['Gene 2', '1.3', '7', '8.12', '1.1', '9.321'],
            ['Gene 3', '6.21', '9.001', '3.76', '6.18', '5']]
    default_observation_lines = [',\t'.join(row) for row in default_observation_matrix]
    # define the default module lines
    default_module_matrix = [
            ['Probe_Set_ID', 'Module', 'Entry Index', 'Average Degree', 'Degree'],
            ['Gene 1', '2', '1', '1.0', '1.0'],
            ['Gene 2', '1', '2', '1.0', '1.0'],
            ['Gene 3', '2', '3', '1.0', '1.0']]
    default_module_lines = [',\t'.join(row) for row in default_module_matrix]
    # define the form objects
    form_objects = [
            Form.MultiLine('observations', 'labeled matrix of observations',
                '\n'.join(default_observation_lines)),
            Form.MultiLine('modules', 'MMC output ',
                '\n'.join(default_module_lines)),
            Form.Float('sigma', 'sigma',
                '0.21', low_exclusive=0),
            Form.RadioGroup('modularity', 'modularity calculation', [
                Form.RadioItem('eric', 'the calculation that eric uses', True),
                Form.RadioItem('other_a', 'my first interpretation'),
                Form.RadioItem('other_b', 'using the formula from the paper'),
                Form.RadioItem('other_c', 'using the formula on wikipedia')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def parse_comma_separated_line(line):
    """
    @param line: a line of comma separated strings
    @return: a list of strings with leading and trailing whitespace removed
    """
    return [value.strip() for value in line.split(',')]

def parse_observation_lines(lines):
    """
    @param lines: lines of comma separated values representing a labeled matrix
    @return: (row labels, column labels, data matrix)
    """
    # get the labeled matrix
    labeled_matrix = [parse_comma_separated_line(line) for line in lines]
    # assert that the matrix is not completely degenerate
    min_nrows = 3
    min_ncols = 3
    if len(labeled_matrix) < min_nrows:
        raise HandlingError('expected at least %d observation rows' % min_nrows)
    if len(labeled_matrix[0]) < min_ncols:
        raise HandlingError('expected at least %d observation columns' % min_ncols)
    # assert that each row has the same number of columns
    ncols = len(labeled_matrix[0])
    for row in labeled_matrix:
        if len(row) != ncols:
            raise HandlingError('expected each observation row to have the same number of columns')
    # get the row labels and the column labels
    row_labels = [row[0] for row in labeled_matrix[1:]]
    column_labels = labeled_matrix[0][1:]
    # get the matrix of data
    data_matrix = []
    for input_row in labeled_matrix[1:]:
        data_row = []
        for element_string in input_row[1:]:
            try:
                element = float(element_string)
            except ValueError, e:
                raise HandlingError('expected each observation element to be a number: ' + element_string)
            data_row.append(element)
        data_matrix.append(data_row)
    return row_labels, column_labels, np.array(data_matrix)

def parse_module_lines(lines):
    """
    @param lines: lines of the MMC csv output
    @return: (gene labels, module indices, gene indices)
    """
    # do some basic validation
    min_nlines = 3
    if len(lines) < min_nlines:
        raise HandlingError('expected at least %d module lines' % min_nlines)
    # extract the parts of the rows of interest
    rows = []
    for line in lines[1:]:
        values = parse_comma_separated_line(line)
        if len(values) != 5:
            raise HandlingError('expected five comma separated values on each module line')
        gene_label, raw_module_index, raw_gene_index, foo, bar = values
        try:
            module_index = int(raw_module_index) - 1
        except ValueError, e:
            raise HandlingError('expected the module index to be an integer: ' + raw_module_index)
        try:
            gene_index = int(raw_gene_index) - 1
        except ValueError, e:
            raise HandlingError('expected the gene index to be an integer: ' + raw_gene_index)
        rows.append([gene_label, module_index, gene_index])
    # return the three lists
    return zip(*rows)

def get_eric_modularity(A, cluster_indices):
    """
    @param A: an affinity matrix defining a graph
    @param cluster_indices: a conformant vector of cluster indices
    """
    # define the number of nodes in the graph and the number of clusters
    n = len(cluster_indices)
    nclusters = max(cluster_indices) + 1
    # initialize some intermediate variables
    volume = 1.0 * sum(sum(A))
    within_cluster = [0] * nclusters
    between_cluster = [0] * nclusters
    # calculate the intermediate variables
    # i and j are node indices
    # a and b are cluster indices
    for i in range(n):
        a = cluster_indices[i]
        for j in range(n):
            b = cluster_indices[j]
            weight = A[i][j]
            between_cluster[a] += weight
            if a == b:
                within_cluster[a] += weight
    # get the modularity from the intermediate variables
    modularity = 0
    for within, between in zip(within_cluster, between_cluster):
        modularity += within/volume - (between/volume)**2
    return modularity

def get_modularity_other_c(A, cluster_indices):
    """
    This is my implementation of modularity using the formula directly from Wikipedia.
    @param A: an affinity matrix defining a graph
    @param cluster_indices: a conformant vector of cluster indices
    """
    # define the number of nodes in the graph and the number of clusters
    n = len(cluster_indices)
    nclusters = max(cluster_indices) + 1
    # define the row sums of the adjacency matrix
    row_sums = [sum(row) for row in A]
    # define one half of the sum of all entries in the adjacency matrix
    m = sum(row_sums) / 2.0
    # define the modularity
    Q = 0
    for i in range(n):
        for j in range(n):
            if cluster_indices[i] == cluster_indices[j]:
                Q += (A[i][j] - row_sums[i] * row_sums[j] / (2*m)) / (2*m)
    return Q

def get_modularity_other_b2(A, cluster_indices):
    """
    This is a modification of the original Girvan-Newman formulation.
    @param A: an affinity matrix defining a graph
    @param cluster_indices: a conformant vector of cluster indices
    """
    # define the number of nodes in the graph and the number of clusters
    n = len(cluster_indices)
    nclusters = max(cluster_indices) + 1
    girvan_e = np.zeros((nclusters, nclusters))
    volume = 0
    for i in range(n):
        for j in range(n):
            if i < j:
                weight = A[i][j]
                volume += weight
                a = cluster_indices[i]
                b = cluster_indices[j]
                if a == b:
                    girvan_e[a][a] += weight
                else:
                    girvan_e[a][b] += weight/2
                    girvan_e[b][a] += weight/2
    for a in range(nclusters):
        for b in range(nclusters):
            girvan_e[a][b] /= volume
    girvan_a = [sum(girvan_e[i]) for i in range(nclusters)]
    modularity = sum(girvan_e[i][i] - girvan_a[i]**2 for i in range(nclusters))
    return modularity

def get_modularity_other_b(A, cluster_indices):
    """
    This is my implementation of modularity using the original Girvan-Newman formulation.
    @param A: an affinity matrix defining a graph
    @param cluster_indices: a conformant vector of cluster indices
    """
    # define the number of nodes in the graph and the number of clusters
    n = len(cluster_indices)
    nclusters = max(cluster_indices) + 1
    girvan_e = np.zeros((nclusters, nclusters))
    volume = 0
    for i in range(n):
        for j in range(n):
            if i < j:
                weight = A[i][j]
                volume += weight
                a = cluster_indices[i]
                b = cluster_indices[j]
                if a == b:
                    girvan_e[a][a] += weight
                else:
                    girvan_e[a][b] += weight
                    girvan_e[b][a] += weight
    for a in range(nclusters):
        for b in range(nclusters):
            girvan_e[a][b] /= volume
    girvan_a = [sum(girvan_e[i]) for i in range(nclusters)]
    modularity = sum(girvan_e[i][i] - girvan_a[i]**2 for i in range(nclusters))
    return modularity

def get_modularity_other_a(A, cluster_indices):
    """
    This was my first implementation of modularity using Eric's definition in his paper.
    @param A: an affinity matrix defining a graph
    @param cluster_indices: a conformant vector of cluster indices
    """
    # define the number of nodes in the graph and the number of clusters
    n = len(cluster_indices)
    nclusters = max(cluster_indices) + 1
    # initialize some intermediate variables
    within_cluster = [0] * nclusters
    between_cluster = [0] * nclusters
    volume = 0
    # calculate the intermediate variables
    # i and j are node indices
    # a and b are cluster indices
    for i in range(n-1):
        a = cluster_indices[i]
        for j in range(i+1, n):
            b = cluster_indices[j]
            weight = A[i][j]
            volume += weight
            if a == b:
                within_cluster[a] += weight
            else:
                between_cluster[a] += weight
                between_cluster[b] += weight
    # get the modularity from the intermediate variables
    modularity = 0
    for within, between in zip(within_cluster, between_cluster):
        modularity += within/volume - ((within+between) / volume)**2
    return modularity

def get_affinity_matrix(corr_matrix, sigma):
    """
    @param corr_matrix: the correlation matrix
    @param sigma: a parameter
    @return: the affinity matrix
    """
    n = len(corr_matrix)
    sigma_squared = sigma * sigma
    affinity_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                # transform correlations to squared distances
                r = corr_matrix[i][j]
                # get the log affinity using sigma and the correlation
                log_affinity = (abs(r) - 1)/sigma_squared
                affinity_matrix[i][j] = math.exp(log_affinity)
    return affinity_matrix

def get_response_content(fs):
    # read the observation lines
    observation_lines = Util.get_stripped_lines(fs.observations.splitlines())
    row_labels, column_labels, data_matrix = parse_observation_lines(observation_lines)
    ngenes = len(row_labels)
    # read the module lines
    module_lines = Util.get_stripped_lines(fs.modules.splitlines())
    gene_labels, module_indices, gene_indices = parse_module_lines(module_lines)
    # each multi-line input should have a header line and N gene lines
    if len(observation_lines) != len(module_lines):
        raise HandlingError('expected the same number of observation lines as MMC output lines')
    # assert that the gene labels match their indices
    for gene_label, gene_index in zip(gene_labels, gene_indices):
        if row_labels[gene_index] != gene_label:
            raise HandlingError('the observation gene order does not appear to match the MMC output order')
    # convert the data matrix to the affinity matrix using sigma
    affinity_matrix = get_affinity_matrix(np.corrcoef(data_matrix), fs.sigma)
    # get the list of module assignments with respect to the observation gene order
    ordered_module_indices = [0] * ngenes
    for module_index, gene_index in zip(module_indices, gene_indices):
        ordered_module_indices[gene_index] = module_index
    # get the modularity according to the specified function
    fn_get_modularity = get_eric_modularity
    if fs.other_a:
        fn_get_modularity = get_modularity_other_a
    elif fs.other_b:
        fn_get_modularity = get_modularity_other_b
    elif fs.other_c:
        fn_get_modularity = get_modularity_other_c
    modularity = fn_get_modularity(affinity_matrix, ordered_module_indices)
    # return the response
    return str(modularity) + '\n'

def main():
    """
    Analyze the Zachary karate club data.
    """
    n = 34
    # create the adjacency matrix
    stripped_lines = Util.get_stripped_lines(g_karate_data.splitlines())
    string_rows = [line.split() for line in stripped_lines if line]
    assert len(string_rows) == n
    for row in string_rows:
        assert len(row) == n
    data_rows = [[float(x) for x in string_row] for string_row in string_rows]
    A = np.array(data_rows)
    # create the ordered module indices
    first_cluster_one_based_indices = [1, 3, 4, 14, 2, 8, 20, 18, 22, 13, 12, 6, 7, 17, 5, 11]
    second_cluster_one_based_indices = [25, 32, 26, 29, 24, 28, 9, 34, 33, 19, 16, 31, 15, 10, 23, 30, 21, 27]
    assert len(first_cluster_one_based_indices + second_cluster_one_based_indices) == n
    assert list(sorted(first_cluster_one_based_indices + second_cluster_one_based_indices)) == range(1, n+1)
    ordered_module_indices = []
    for i in range(n):
        if i+1 in first_cluster_one_based_indices:
            ordered_module_indices.append(0)
        else:
            ordered_module_indices.append(1)
    # print the modularity
    Q = get_modularity_other_b(A, ordered_module_indices)
    print 'modularity calculated using my interpretation of the method of the paper', Q
    Q = get_modularity_other_b2(A, ordered_module_indices)
    print 'modularity calculated using a modification of my interpretation of the method of the paper', Q
    Q = get_modularity_other_c(A, ordered_module_indices)
    print 'modularity calculated using the method on wikipedia', Q
    Q = get_eric_modularity(A, ordered_module_indices)
    print 'modularity calculated using the method eric used:', Q
    print 'expected modularity: .375 +/- .025'

if __name__ == '__main__':
    main()
