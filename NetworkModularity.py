"""
Implement the Newman and Girvan calculation of network modularity.
This objective function defines how well a given partition of a graph
represents its underlying connectivity.
"""

def get_eric_modularity(A, cluster_indices):
    """
    @param A: an affinity matrix defining a graph
    @param cluster_indices: a conformant vector of cluster indices
    """
    # define the number of nodes in the graph and the number of clusters
    n = len(cluster_indices)
    nclusters = max(cluster_indices) + 1
    # initialize some intermediate variables
    volume = 1.0 * sum(sum(row) for row in A)
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
    girvan_e = numpy.zeros((nclusters, nclusters))
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
    girvan_e = numpy.zeros((nclusters, nclusters))
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

def main():
    A = [
            [0, 1, 1, 0, 0, 0],
            [1, 0, 1, 0, 0, 0],
            [1, 1, 0, 1, 0, 0],
            [0, 0, 1, 0, 1, 1],
            [0, 0, 0, 1, 0, 1],
            [0, 0, 0, 1, 1, 0]]
    cluster = [0, 0, 0, 1, 1, 1]
    print 'testing the symmetric example'
    Q = get_eric_modularity(A, cluster)
    print 'calculated using the method of eric:', Q
    print 'calculated by hand:', 5.0 / 14.0

if __name__ == '__main__':
    main()
