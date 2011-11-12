"""
Given N free energies, get an NxN rate matrix.

Updated 20111111.
"""

from StringIO import StringIO
import math
import numpy as np
import scipy

import Form
import FormOut
import iterutils
from MatrixUtil import ndot

#TODO maybe this should return a rate matrix object instead of a report

def get_form():
    """
    @return: the body of a form
    """
    # define the default energy string
    default_energies = [2, 4, 6, 8]
    default_energy_string = '\n'.join(str(E) for E in default_energies)
    # define the form objects
    form_objects = [
            Form.MultiLine('energies', 'ordered energies',
                default_energy_string)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def R_to_distn(R):
    """
    @param R: rate matrix
    @return: stationary distribution
    """
    n = len(R)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    dummy, pi_eigenvector = min(val_vec_pairs)
    total = np.sum(pi_eigenvector)
    pi_arr = np.array([v/total for v in pi_eigenvector])
    return pi_arr

def get_numerical_identicality(R, t):
    """
    @param R: rate matrix
    @param t: time
    @return: probability that two t-separated observations are identical
    """
    pi_arr = R_to_distn(R)
    T = scipy.linalg.expm(R*t)
    return np.dot(pi_arr, np.diag(T))

def get_identicality_params(R):
    """
    This returns the parameters for an identicality function.
    If the rate matrix has n states
    then the identicality function is
    f(t) = a1*exp(b1*t) + a2*exp(b2*t) + ... + a{n-1}*exp(b{n-1}*t) + c
    @param R: time reversible rate matrix
    @return: a array, b array, c
    """
    n = len(R)
    pi_arr = R_to_distn(R)
    # symmetrize
    lam = np.diag(np.sqrt(pi_arr))
    rlam = np.diag(np.reciprocal(np.sqrt(pi_arr)))
    S = ndot(lam, R, rlam)
    print 'S should be symmetric:'
    print S
    print S - S.T
    # eigendecompose the symmetric matrix
    W, V = scipy.linalg.eigh(S)
    w_v_pairs = [(W[i], V[:,i]) for i in range(n)]
    # get the exponential coefficients
    eps = 1e-12
    identicality_coeffs = [
            np.dot(pi_arr, v*v) for w, v in w_v_pairs if abs(w) > eps]
    # get the exponential rate constants
    identicality_rates = [
            w for w in W if abs(w) > eps]
    # get the one dimensional constant
    identicality_const = np.inner(pi_arr, pi_arr)
    # return the identicality parameters
    return (identicality_coeffs, identicality_rates, identicality_const)

def get_symbolic_identicality(coeffs, rates, c, t):
    """
    @param coeffs: exponential coefficients
    @param rates: exponential rate constants
    @param c: a numerical constant
    @param t: elapsed time
    """
    return c + sum(coeff * math.exp(r*t) for coeff, r in zip(coeffs, rates))

def get_identicality_derivative(coeffs, rates, t):
    return sum(r * coeff * math.exp(r*t) for coeff, r in zip(coeffs, rates))


def get_response_content(fs):
    # read the energies from the form data
    energies = []
    for line in iterutils.stripped_lines(fs.energies.splitlines()):
        try:
            energy = float(line)
        except ValueError as e:
            raise ValueError('invalid energy: %s' % line)
        energies.append(energy)
    n = len(energies)
    if n > 100:
        raise ValueError('too many energies')
    # compute the rate matrix
    R = np.zeros((n, n))
    for row in range(n):
        for col in range(n):
            rate = math.exp(-(energies[col] - energies[row]))
            R[row, col] = rate
    for i, r in enumerate(R):
        R[i, i] = -np.sum(r) + 1
    # get the transition matrix at large finite time
    large_t = 1000.0
    T = scipy.linalg.expm(R*large_t)
    # eigendecompose
    Wr, Vr = scipy.linalg.eig(R, left=False, right=True)
    Wl, Vl = scipy.linalg.eig(R, left=True, right=False)
    # get left eigenvector associated with stationary distribution
    val_vec_pairs = [(abs(Wl[i]), Vl[:,i]) for i in range(n)]
    dummy, pi_eigenvector = min(val_vec_pairs)
    # get the stationary distribution itself
    total = np.sum(pi_eigenvector)
    pi_arr = np.array([v/total for v in pi_eigenvector])
    # get the square root stationary vector and diagonal matrix
    sqrt_pi_arr = np.sqrt(pi_arr)
    lam = np.diag(sqrt_pi_arr)
    # get reciprocal arrays
    recip_sqrt_pi_arr = np.reciprocal(sqrt_pi_arr)
    recip_lam = np.reciprocal(lam)
    # print things
    np.set_printoptions(linewidth=300)
    out = StringIO()
    print >> out, 'rate matrix:'
    print >> out, R
    print >> out
    print >> out, 'rate matrix row sums:'
    print >> out, np.sum(R, axis=1)
    print >> out
    print >> out, 'eigenvalues:'
    print >> out, Wr
    print >> out
    print >> out, 'corresponding orthonormal right eigenvectors (columns):'
    print >> out, Vr
    print >> out
    print >> out, 'eigenvalues:'
    print >> out, Wl
    print >> out
    print >> out, 'corresponding orthonormal left eigenvectors (columns):'
    print >> out, Vl
    print >> out
    print >> out, 'L2 normalized eigenvector associated with stationary distn:'
    print >> out, pi_eigenvector
    print >> out
    print >> out, 'L1 renormalized vector (the stationary distribution):'
    print >> out, pi_arr
    print >> out
    print >> out
    # eigendecompose the transition matrix
    Wr, Vr = scipy.linalg.eig(T, left=False, right=True)
    Wl, Vl = scipy.linalg.eig(T, left=True, right=False)
    print >> out, 'transition matrix for t=%f:' % large_t
    print >> out, T
    print >> out
    print >> out, 'transition matrix row sums:'
    print >> out, np.sum(T, axis=1)
    print >> out
    print >> out, 'eigenvalues:'
    print >> out, Wr
    print >> out
    print >> out, 'corresponding orthonormal right eigenvectors (columns):'
    print >> out, Vr
    print >> out
    print >> out, 'eigenvalues:'
    print >> out, Wl
    print >> out
    print >> out, 'corresponding orthonormal left eigenvectors (columns):'
    print >> out, Vl
    print >> out
    print >> out, 'incorrect reconstitution of the transition matrix:'
    print >> out, ndot(Vr, np.diag(Wr), Vl.T)
    print >> out
    print >> out
    # Use the known properties of reversibility to symmetrize the matrix.
    t = 3
    coeffs, rates, c = get_identicality_params(R)
    print >> out, 'brute identicality computation for t=%f:' % t
    print >> out, get_numerical_identicality(R, t)
    print >> out
    print >> out, 'sophisticated identicality computation for t=%f:' % t
    print >> out, get_symbolic_identicality(coeffs, rates, c, t)
    print >> out
    print >> out
    # Try another couple rate matrices.
    e2 = math.exp(2)
    en2 = math.exp(-2)
    rate_matrices = [
            np.array([[-2.0, 2.0], [2.0, -2.0]]),
            np.array([[-1.0, 1.0], [3.0, -3.0]]),
            np.array([[-1, 1, 0], [1, -2, 1], [0, 1, -1]]),
            #np.array([[-4.0, 4.0, 0], [1, -2, 1], [0, 4, -4]])]
            #np.array([[-1, 1, 0], [7, -14, 7], [0, 1, -1]])]
            np.array([[-en2, en2, 0], [e2, -2*e2, e2], [0, en2, -en2]])]
    t = 3.0
    for R in rate_matrices:
        coeffs, rates, c = get_identicality_params(R)
        print >> out, 'test rate matrix:'
        print >> out, R
        print >> out
        print >> out, 'eigenvalues:'
        print >> out, scipy.linalg.eigvals(R)
        print >> out
        print >> out, 'stationary distribution:'
        print >> out, R_to_distn(R)
        print >> out
        print >> out, 'brute identicality computation for t=%f:' % t
        print >> out, get_numerical_identicality(R, t)
        print >> out
        print >> out, 'sophisticated identicality computation for t=%f:' % t
        print >> out, get_symbolic_identicality(coeffs, rates, c, t)
        print >> out
        print >> out, 'identicality derivative for t=%f:' % t
        print >> out, get_identicality_derivative(coeffs, rates, t)
        print >> out
        print >> out
    # return the message
    return out.getvalue().rstrip()

