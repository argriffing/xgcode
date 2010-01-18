"""
Check spectral splits of some small perturbations
distance matrices of symmetric trees.
"""

import numpy as np

import Euclid

def show_split(D):
    HSH = Euclid.edm_to_dccov(D)
    # get the eigendecomposition
    eigenvalues, V_T = np.linalg.eigh(HSH)
    eigenvectors = V_T.T.tolist()
    # get the eigenvalue and eigenvector of interest
    w, v = max(zip(eigenvalues, eigenvectors))
    # show the results
    print 'the maximum of these eigenvalues is interesting:'
    print '\t'.join(str(x) for x in sorted(eigenvalues))
    print 'the interesting eigenvector:'
    print '\t'.join(str(x) for x in v)

def do_a():
    D = np.array([
        [0, 2, 8, 8, 8, 8],
        [2, 0, 8, 8, 8, 8],
        [8, 8, 0, 2, 8, 8],
        [8, 8, 2, 0, 8, 8],
        [8, 8, 8, 8, 0, 2],
        [8, 8, 8, 8, 2, 0]], dtype=np.float64)
    print 'trisymmetric results:'
    show_split(D)
    print

def do_b():
    D = np.array([
        [0, 2, 7, 7, 7, 7],
        [2, 0, 7, 7, 7, 7],
        [7, 7, 0, 2, 8, 8],
        [7, 7, 2, 0, 8, 8],
        [7, 7, 8, 8, 0, 2],
        [7, 7, 8, 8, 2, 0]], dtype=np.float64)
    print 'bisymmetric results:'
    show_split(D)
    print

def do_perturbed_bisymmetric_matrix(perturbation):
    d = perturbation
    D = np.array([
        [0, 2+d, 7-d, 7-d, 7+d, 7+d],
        [2+d, 0, 7+d, 7+d, 7-d, 7-d],
        [7-d, 7+d, 0, 2-d, 8+d, 8+d],
        [7-d, 7+d, 2-d, 0, 8+d, 8+d],
        [7+d, 7-d, 8+d, 8+d, 0, 2-d],
        [7+d, 7-d, 8+d, 8+d, 2-d, 0]], dtype=np.float64)
    print 'bisymmetric results perturbed by', d, ':'
    show_split(D)
    show_analytical_perturbed_split(perturbation)
    print

def show_analytical_perturbed_split(perturbation):
    # get the analytical values
    eigenvalues = list(gen_analytical_eigenvalues(perturbation))
    v = normalized(gen_analytical_eigenvector(perturbation))
    # show the results
    print 'analytical eigenvalues:'
    print '\t'.join(str(x) for x in sorted(eigenvalues))
    print 'analytical eigenvector:'
    print '\t'.join(str(x) for x in v)

def gen_analytical_eigenvalues(perturbation):
    d = perturbation
    a = (9*d*d + 12*d + 36)**.5
    yield -a/2 + d + 4
    yield a/2 + d + 4
    yield 17/3. - d/2.
    yield 1 - d/2.
    yield 1 - d/2.
    yield 0.

def gen_analytical_eigenvector(perturbation):
    d = perturbation
    a = ((9*d*d + 12*d + 36)**.5 + d + 6 ) / (4*d)
    yield 1.
    yield -1.
    yield -a
    yield -a
    yield a
    yield a

def normalized(v_in):
    v = list(v_in)
    mag = sum(x**2 for x in v)**.5
    return [x / mag for x in v]

def main():
    do_b()
    for i in range(4):
        do_perturbed_bisymmetric_matrix(10**-i)

if __name__ == '__main__':
    main()

