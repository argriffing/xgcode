"""
Use the standard k-means clustering algorithm.

Use an example from the clusterSim R package.
"""

import unittest
import random
import textwrap
import collections

import numpy as np

import iterutils
import MatrixUtil
import Util
import Form

g_data_ratio = textwrap.dedent("""
         v_1       v_2       v_3       v_4       v_5
         1   4.790679  5.078952  8.485110 10.430934 13.293767
         2   5.045025  4.951188 11.628062  9.599697 11.739800
         3   4.729506  5.573883  8.664898  9.898832 15.153742
         4   5.106198  4.456256 11.448274 10.131799  9.879825
         5   6.349961  5.233782  5.437522  8.766744 12.901515
         6   3.485743  4.796357 14.675650 11.263887 12.132052
         7   4.767566  3.620416  8.896164 12.082252  7.801038
         8   5.068138  6.409723 11.217008  7.948379 17.232529
         9   5.100619  5.993950 10.490193  8.537678 15.836951
         10  4.735085  4.036190  9.622979 11.492953  9.196616
         11  5.799416  6.789289 11.706279  6.420144 17.702119
         12  4.036288  3.240851  8.406893 13.610487  7.331448
         13  6.881446  3.985094  7.253590  9.140405  7.213446
         14  2.954258  6.045045 12.859582 10.890227 17.820121
         15  6.150072  4.789024  9.356354  8.760574 10.480844
         16  5.244338 13.564942 10.756818 11.270057 14.552723
         17  4.719234 16.503408 13.959299  8.021049 14.129663
         18  4.557217 14.505816  6.153873 12.009583 10.903904
         19  5.406355 15.562534  9.815691 12.701430 16.534460
         20  4.692166 15.322574 10.297480  7.329201  8.499107
         21  5.271406 14.745776 10.936124  9.189059 13.491379
         22  4.375457 13.267942  9.177048 10.841572 11.542188
         23  5.588115 16.800408 12.344995 10.167320  6.546140
         24  2.539172 15.090441  7.768177  9.863311 18.487427
         25  7.424400 14.977909 10.329711 13.081164 14.476322
         26  4.837496 14.044788  9.783460  6.949467 10.557245
         27  5.126075 16.023562  8.938004 10.033222 12.645054
         28  5.863418 14.165410 11.175168  9.997409 12.388513
         29  4.100154 15.902940  9.209020  5.659804 13.320272
         30  5.254999 14.521551 10.904152 14.370827 11.713295
         31  9.952174 10.665441 11.263631  9.329913 11.941795
         32  9.996607  9.160925  8.849541 10.700718 13.091772
         33 10.095053 10.626930 15.285856  9.817046 14.052205
         34  9.853728  9.199436  4.827316 10.213585 10.981362
         35  9.291071 11.390109 13.334334 12.713901 16.902076
         36 10.657710  8.436257  6.778838  7.316730  8.131491
         37  8.588608 11.231346  8.982902  9.867728 12.583535
         38 11.360173  8.595020 11.130270 10.162903 12.450032
         39  9.520658  9.087109 11.760918  8.505001 12.893298
         40 10.428123 10.739257  8.352254 11.525630 12.140269
         41  7.991082  8.915160  8.072875 12.002505 13.603081
         42 11.957699 10.911206 12.040297  8.028126 11.430486
         43  9.031321  9.723553  5.954012 13.282765 14.357532
         44 10.917460 10.102813 14.159160  6.747866 10.676035
         45 10.358533 11.215437  9.343065  5.029610 13.610742
         46 16.526776  5.993841 10.770106 15.001021 11.422825
         47 13.251358  4.132012  7.553808 11.825198  9.348163
         48 14.912288  6.259952 12.559364  8.205433 15.685404
         49 14.865847  3.865901 10.411889 13.458668 14.936874
         50 15.551555  4.639883  9.701283  6.571963 10.096693
         51 14.226580  5.485970  8.590946 12.142722 14.003842
         52 13.235911  5.132797 11.522226  7.887909 11.029725
         53 16.542223  4.993056  7.295439 13.491249 14.359622
         54 14.872899  5.281239 12.817733  6.539382 10.673945
         55 14.905235  4.844615 14.415468  6.948218 13.263033
         56 14.765832  7.179090  5.697703 13.082413 11.770534
         57 15.012302  2.946763 10.701978  9.923576 10.465631
         58 15.470419  4.788885  9.411194 10.107055 14.567936
         59 14.307716  5.336968 12.564805 12.304334 13.132195
         60 16.553059  4.119027  7.548367  7.726297 11.901372
         61 14.044068 16.425270 11.893600 10.249055  7.908837
         62 16.052005 13.567111  8.219572  9.781576 17.124730
         63 14.994919 13.249362  5.228936  7.461541 12.564844
         64 15.101153 16.743019 14.884236 12.569090 12.468723
         65 16.655148 15.142145 10.942706  8.800695 12.669968
         66 13.440925 14.850236  9.170466 11.229936 12.363599
         67 14.222441 15.832170  7.002342 10.938726 12.113179
         68 15.873632 14.160211 13.110829  9.091905 12.920388
         69 15.226869 15.846884  7.310058  7.110427  7.016356
         70 14.869203 14.145497 12.803114 12.920205 18.017211
         71 14.766368 14.691712  2.181596 13.562787  9.280304
         72 15.329705 15.300669 17.931576  6.467844 15.753263
         73 16.742500 15.927990 12.023954  7.113512 16.367873
         74 13.353573 14.064391  8.089218 12.917119  8.665694
         75 14.327491 15.053332  5.812640  8.866648 11.258022
         """).strip()

g_cl = textwrap.dedent("""
   1  1  2  1  2  1  3  2  2  3  2  3  3
   2  1  4  4  4  4  4  4  4  4  4  4  4 
   4  4  4  4  5  6  5  6  5  6  6  5  5
   6  6  5  6  5  5  7  7  8  7  8  7  8 
   7  8  8  7  8  7  7  8  9  9 10  9  9
   10 10  9 10  9 10  9  9 10 10
   """).strip()

g_expected_calinski = 25.15651


class InitStrategy(Form.RadioGroup):
    
    def __init__(self):
        Form.RadioGroup.__init__(
                self, 'kmeans_init', 'cluster center initialization', [
                    Form.RadioItem('init_choice',
                        'center on randomly chosen observed points', True),
                    Form.RadioItem('init_cluster',
                        'use centroids of a random partition'),
                    Form.RadioItem('init_range',
                        'choose centers uniformly in the observed range')])

    def add_argument(self, parser):
        """
        Add an argument to the argparse parser.
        """
        return parser.add_argument('--kmeans_init', default='init_choice',
                choices=['init_choice', 'init_cluster', 'init_range'],
                help='cluster center initialization')

    def string_to_function(self, kmeans_init):
        d = {
                'init_cluster' : gen_random_centers_via_clusters,
                'init_choice' : gen_random_centers_via_choice,
                'init_range' : gen_random_centers_via_range}
        return d[kmeans_init]


def gen_random_centers_via_range(points, nclusters):
    """
    Yield guesses for initial clusters.
    This methods is not so great.
    """
    mymin = np.min(points, axis=0)
    mymax = np.max(points, axis=0)
    for i in range(nclusters):
        yield np.array([random.uniform(a, b) for a, b in zip(mymin, mymax)])

def gen_random_centers_via_choice(points, nclusters):
    """
    Return a list of guesses for initial clusters.
    This is the preferred way to get initial guesses.
    """
    return random.sample(points, nclusters)

def gen_random_centers_via_clusters(points, nclusters):
    """
    Return a list of guesses for initial clusters.
    This methods is not so great.
    """
    # get random labels with each label appearing at least once
    npoints = len(points)
    labels = np.random.randint(0, nclusters, npoints)
    indices = random.sample(range(npoints), nclusters)
    for index, label in zip(indices, range(nclusters)):
        labels[index] = label
    return get_centers(points, labels)

def get_point_center_sqdists(points, centers):
    """
    Inputs and outputs are numpy arrays.
    @param points: the points to be reclustered
    @param centers: cluster centers
    @return: for each point, the squared distance to each center
    """
    MatrixUtil.assert_2d(points)
    MatrixUtil.assert_2d(centers)
    npoints = len(points)
    ncenters = len(centers)
    # get the dot products of points with themselves
    pself = np.sum(points*points, axis=1)
    # get the dot products of centers with themselves
    cself = np.sum(centers*centers, axis=1)
    # get the matrix product of points and centers
    prod = np.dot(points, centers.T)
    # get the matrix of squared distances
    sqdists = (
        np.outer(pself, np.ones(ncenters)) +
        np.outer(np.ones(npoints), cself) -
        2*prod)
    return sqdists

def get_centers(points, labels):
    """
    Inputs and outputs are numpy arrays.
    @param points: euclidean points
    @param labels: conformant cluster indices
    """
    MatrixUtil.assert_2d(points)
    MatrixUtil.assert_1d(labels)
    if len(points) != len(labels):
        raise ValueError('array incompatibility')
    ncoords = len(points[0])
    nclusters = max(labels) + 1
    sums = [np.zeros(ncoords) for i in range(nclusters)]
    counts = [0]*nclusters
    for point, label in zip(points, labels):
        sums[label] += point
        counts[label] += 1
    M = np.array([s/c for s, c in zip(sums, counts)])
    return M

def get_labels_without_cluster_removal(sqdists):
    """
    Inputs and outputs are numpy arrays.
    @param sqdists: for each point, the squared distance to each center
    @return: for each point, the label of the nearest cluster
    """
    return np.argmin(sqdists, axis=1)

def get_labels(sqdists):
    """
    Inputs and outputs are numpy arrays.
    Account for the fact that sometimes a cluster will go away.
    That is, if no point is in the voronoi region of a centroid,
    then in the next iteration this cluster should disappear.
    @param sqdists: for each point, the squared distance to each center
    @return: for each point, the label of the nearest cluster
    """
    labels = np.argmin(sqdists, axis=1)
    new_to_old = list(iterutils.unique_everseen(labels))
    old_to_new = dict((old, new) for new, old in enumerate(new_to_old))
    return np.array([old_to_new[old] for old in labels])

def get_wcss(sqdists, labels):
    """
    Get the within-cluster sum of squares.
    @param sqdists: for each point, the squared distance to each center
    @param labels: cluster labels
    @return: within-cluster sum of squares
    """
    MatrixUtil.assert_2d(sqdists)
    MatrixUtil.assert_1d(labels)
    if len(sqdists) != len(labels):
        raise ValueError('array incompatibility')
    return sum(row[label] for row, label in zip(sqdists, labels))

def lloyd(points, labels):
    """
    This is the standard algorithm for kmeans clustering.
    @param points: points in euclidean space
    @param labels: initial cluster labels
    @return: within cluster sum of squares, and labels
    """
    while True:
        centers = get_centers(points, labels)
        sqdists = get_point_center_sqdists(points, centers)
        next_labels = get_labels(sqdists)
        if np.array_equal(next_labels, labels):
            wcss = get_wcss(sqdists, labels)
            return wcss, labels
        labels = next_labels

def lloyd_with_restarts(points, nclusters, nrestarts, init_strategy):
    """
    This is the standard algorithm for kmeans clustering with restarts.
    @param points: points in euclidean space
    @param nclusters: the number of clusters
    @param nrestarts: the number of random restarts
    @param init_strategy: a function that guesses initial clusters
    @return: within cluster sum of squares, labels
    """
    npoints = len(points)
    best_wcss = None
    best_labels = None
    for i in range(nrestarts):
        centers = np.array(list(
            init_strategy(points, nclusters)))
        sqdists = get_point_center_sqdists(points, centers)
        labels = get_labels(sqdists)
        wcss, labels = lloyd(points, labels)
        if (best_wcss is None) or (wcss < best_wcss):
            best_wcss = wcss
            best_labels = labels
    return best_wcss, best_labels

def get_allmeandist(points):
    """
    Use this to find the bgss argument for the calinski index.
    @param points: each row is a point and each column is a coordinate
    @return: a single number
    """
    return np.sum((points - np.mean(points, axis=0))**2)

def get_calinski_index(bgss, wgss, k, n):
    """
    Compute the calinski index.
    This is a criterion that can be used to pick the number of clusters;
    the number of clusters that gives the greatest index could be used.
    @param bgss: between groups sum of squares
    @param wgss: within groups sum of squares
    @param k: number of clusters
    @param n: number of points
    @return: a floating point number
    """
    if not (1 < k < n):
        raise ValueError(
                'the calinski index '
                'is defined for integers k and n such that 1 < k < n')
    numerator = bgss / float(k - 1)
    denominator = wgss / float(n - k)
    if not denominator:
        return float('inf')
    return numerator / denominator

def get_calinski_index_naive(points, labels):
    """
    Compute the calinski index.
    This is a naive computation for testing.
    @param points: points in euclidean space
    @param labels: a conformant array giving a cluster label to each point
    """
    # convert the labels to a cluster map
    cluster_map = collections.defaultdict(set)
    for i, label in enumerate(labels):
        cluster_map[label].add(i)
    # get the sum of squared distances to the center
    allmeandist = get_allmeandist(points)
    # check the input sizes
    n = len(points)
    k = len(cluster_map)
    if not (1 < k < n):
        raise ValueError(
                'the calinski index '
                'is defined for integers k and n such that 1 < k < n')
    # compute wgss using squares of distances
    wgss = 0
    for point_indices in cluster_map.values():
        d = 0
        for i in point_indices:
            for j in point_indices:
                diff = points[j] - points[i]
                d += np.dot(diff, diff)
        wgss += d / float(2*len(point_indices))
    # get the index from allmeandist, wgss, and the sizes
    bgss = allmeandist - wgss
    numerator = bgss / float(k - 1)
    denominator = wgss / float(n - k)
    if not denominator:
        return float('inf')
    return numerator / denominator


class TestCalinski(unittest.TestCase):

    def create_test_points(self):
        data_lines = Util.get_stripped_lines(g_data_ratio.splitlines())[1:]
        data_rows = [x.split()[1:] for x in data_lines]
        points = [[float(x) for x in row] for row in data_rows]
        return np.array(points)

    def create_test_labels(self):
        line = ' '.join(g_cl.splitlines())
        row = line.split()
        labels = [int(x) for x in row]
        return np.array(labels)

    def test_wcss_bcss_regression(self):
        points = np.array([[0],[1],[4]], dtype=float)
        labels = np.array([0,0,1])
        centers = get_centers(points, labels)
        sqdists = get_point_center_sqdists(points, centers)
        allmeandist = get_allmeandist(points)
        wcss = get_wcss(sqdists, labels)
        bcss = allmeandist - wcss
        self.assertTrue(np.allclose(wcss, 0.5))
        self.assertTrue(np.allclose(bcss, 49 / 6.0))

    def test_calinski(self):
        points = self.create_test_points()
        labels = self.create_test_labels()
        centers = get_centers(points, labels)
        sqdists = get_point_center_sqdists(points, centers)
        allmeandist = get_allmeandist(points)
        wcss = get_wcss(sqdists, labels)
        bcss = allmeandist - wcss
        k = len(set(labels))
        n = len(points)
        calinski = get_calinski_index(bcss, wcss, k, n)
        fcal_observed = '%.4f' % calinski
        fcal_expected = '%.4f' % g_expected_calinski
        self.assertEqual(fcal_observed, fcal_expected)

    def test_calinski_naive(self):
        points = self.create_test_points()
        labels = self.create_test_labels()
        calinski = get_calinski_index_naive(points, labels)
        fcal_observed = '%.4f' % calinski
        fcal_expected = '%.4f' % g_expected_calinski
        self.assertEqual(fcal_observed, fcal_expected)


if __name__ == '__main__':
    unittest.main()
