"""
Attempt to find branch lengths to cause a 2D MDS collision.

All vertices must have integer labels.
"""

from StringIO import StringIO
import math
import time

import numpy as np
import scipy
from scipy import linalg
from scipy import optimize

import Form
import FormOut
import Ftree
import FtreeIO


def get_form():
    """
    @return: the body of a form
    """
    # define default tree strings
    true_s = '((1:3, 2:6)5:1, 3:9, 4:12)6;'
    test_s = '((1, 3)5, 2, 4)6;'
    # define the form objects
    form_objects = [
            Form.MultiLine('true_tree', 'true tree', true_s),
            Form.MultiLine('test_tree', 'test topology', test_s)]
    return form_objects

def get_form_out():
    return FormOut.Report()


class Functor:

    def __init__(self, T_test, Vp, C, w):
        """
        @param T_test: test topology
        @param Vp: matrix of column eigenvectors of Schur complement in L
        @param C: target matrix
        @param w: Fiedler and Fiedler+1 eigenvalues Schur complement in L
        """
        self.T_test = T_test
        self.Vp = Vp
        self.C = C
        self.w = w
        # define an edge order
        self.u_edges = sorted(self.T_test)
        # precompute vertex lists
        self.leaves = Ftree.T_to_leaves(T_test)
        self.internal = Ftree.T_to_internal_vertices(T_test)
        self.vertices = self.leaves + self.internal

    def X_to_B_Vr(self, X):
        """
        Unpack purely from X.
        """
        start = 0
        # extract log branch lengths
        n = len(self.u_edges)
        log_branch_lengths = X[start:start+n]
        start += n
        # extract first leaf eigenvector
        n = len(self.internal)
        vr1 = X[start:start+n]
        start += n
        # extract second leaf eigenvector
        n = len(self.internal)
        vr2 = X[start:start+n]
        start += n
        # build some more complicated data structures
        B = dict((u_edge, math.exp(logb)) for u_edge, logb in zip(
            self.u_edges, log_branch_lengths))
        Vr = np.vstack([vr1, vr2]).T
        # return the unpacked values
        return B, Vr

    def X_to_L_V(self, X):
        """
        Unpack in a way that uses initialized state.
        """
        B, Vr = self.X_to_B_Vr(X)
        # get the laplacian matrix
        L = Ftree.TB_to_L_principal(self.T_test, B, self.vertices)
        # get the augmented vector
        V = np.vstack([self.Vp, Vr])
        # return the unpacked values
        return L, V

    def __call__(self, X):
        """
        First few entries of X are logs of branch lengths.
        Next few entries are vr1 entries.
        Next few entries are vr2 entries.
        @param X: a 1D numpy array of floats
        """
        # unpack the parameter array
        B, Vr = self.X_to_B_Vr(X)
        L, V = self.X_to_L_V(X)
        # get the error matrix
        E = np.dot(L, V) - self.C
        # compute the squared frobenius norm of the error
        frob_norm_err = np.sum(E*E)
        # use a hack to make sure we are really using the first two eigenvalues
        L_schur = Ftree.TB_to_L_schur(self.T_test, B, self.leaves)
        w_observed = scipy.linalg.eigvalsh(L_schur, eigvals=(1,2))
        w_error = w_observed - self.w
        eigenvalue_err = np.sum(w_error*w_error)
        # return the total error
        return frob_norm_err + eigenvalue_err


def get_response_content(fs):
    nseconds_limit = 5.0
    R_true, B_true = FtreeIO.newick_to_RB(fs.true_tree, int)
    R_test = FtreeIO.newick_to_R(fs.test_tree, int)
    # get the unrooted tree topology
    T_true = Ftree.R_to_T(R_true)
    T_test = Ftree.R_to_T(R_test)
    # check the trees for vertex compatibility
    if set(Ftree.T_to_order(T_true)) != set(Ftree.T_to_order(T_test)):
        raise ValueError('vertex sets are not equal')
    if set(Ftree.T_to_leaves(T_true)) != set(Ftree.T_to_leaves(T_test)):
        raise ValueError('leaf vertex sets are not equal')
    if set(Ftree.T_to_internal_vertices(T_true)) != set(
            Ftree.T_to_internal_vertices(T_test)):
        raise ValueError('internal vertex sets are not equal')
    # get the 2D MDS for the true tree
    leaves = Ftree.T_to_leaves(T_true)
    internal = Ftree.T_to_internal_vertices(T_true)
    vertices = leaves + internal
    L_schur = Ftree.TB_to_L_schur(T_true, B_true, leaves)
    w_all, Vp_all = scipy.linalg.eigh(L_schur)
    w, Vp = w_all[1:3], Vp_all[:, 1:3]
    # make the constant matrix for Frobenius norm comparison
    C = np.zeros((len(vertices), 2))
    C[:len(leaves)] = w*Vp
    # keep doing iterations until we run out of time
    mymax = 256
    t_initial = time.time()
    while time.time() - t_initial < nseconds_limit / 2:
        mymax *= 2
        f = Functor(T_test.copy(), Vp.copy(), C.copy(), w.copy())
        initial_guess = np.ones(len(T_test) + 2*len(internal))
        results = scipy.optimize.fmin(
                f, initial_guess, ftol=1e-8, xtol=1e-8, full_output=True,
                maxfun=mymax, maxiter=mymax)
        xopt, fopt, itr, funcalls, warnflag = results
    # look at the values from the longest running iteration
    B, Vr = f.X_to_B_Vr(xopt)
    L, V = f.X_to_L_V(xopt)
    Lrr = Ftree.TB_to_L_block(T_test, B, internal, internal)
    Lrp = Ftree.TB_to_L_block(T_test, B, internal, leaves)
    H_ext = -np.dot(np.linalg.pinv(Lrr), Lrp)
    N = dict((v, str(v)) for v in vertices)
    # start writing the response
    out = StringIO()
    print >> out, 'xopt:', xopt
    print >> out, 'fopt:', fopt
    print >> out, 'number of iterations:', itr
    print >> out, 'number of function calls:', funcalls
    print >> out, 'warning flags:', warnflag
    print >> out, 'first four eigenvalues:', w_all[:4]
    print >> out, 'Vr:'
    print >> out, Vr
    print >> out, '-Lrr^-1 Lrp Vp:'
    print >> out, np.dot(H_ext, Vp)
    print >> out, C
    print >> out, np.dot(L, V)
    print >> out, FtreeIO.RBN_to_newick(R_test, B, N)
    return out.getvalue()

