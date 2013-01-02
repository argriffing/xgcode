"""
Guess rational equilibrium distributions associated with Moran processes.
"""

from StringIO import StringIO
import itertools

import numpy as np
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

import Form
import FormOut
from MatrixUtil import ndot

def get_form():
    """
    @return: the body of a form
    """
    return [
            #Form.Integer('N', 'population size', 3, low=1, high=5),
            ]

def get_form_out():
    return FormOut.Report()



##############################################################################
# These functions define an ordering on vertices in a lattice in a tetrahedron.

def gen_states(N, k):
    """
    Yield multinomial states.
    Each state is a list of length k and with sum N.
    The state ordering is one which simple to generate recursively.
    @param N: population size
    @param k: number of bins
    """
    if k == 1:
        yield [N]
    elif k == 2:
        for i in range(N+1):
            yield [i, N-i]
    else:
        for i in range(N+1):
            for suffix in gen_states(N-i, k-1):
                yield [i] + suffix

def get_inverse_dict(M):
    """
    The input M[i,j] is count of allele j for pop state index i.
    The output T[(i,j,...)] maps allele count tuple to pop state index
    @param M: multinomial state map
    @return: T
    """
    return dict((tuple(state), i) for i, state in enumerate(M))


##############################################################################
# These functions construct a sparse rate matrix.

def _coo_helper(coo_i, coo_j, coo_r, i, j, r):
    """
    Update the lists which will be used to construct a scipy.sparse.coo_matrix.
    @param coo_i: list of source state indices
    @param coo_j: list of sink state indices
    @param coo_r: list of rates
    @param i: source state index
    @param j: sink state index
    @param r: rate from source to sink
    """
    # add to rate from source to sink
    coo_i.append(i)
    coo_j.append(j)
    coo_r.append(r)
    # add to rate out of the source
    coo_i.append(i)
    coo_j.append(i)
    coo_r.append(-r)

def create_coo_moran(M, T, alpha):
    """
    Construct the sparse Moran rate matrix.
    Initially build it as a coo_matrix,
    which presumably is efficient to transpose and to convert to a csc_matrix.
    @param M: index to state description
    @param T: state description to index
    @return: scipy.sparse.coo_matrix
    """
    nstates = M.shape[0]
    ci = []
    cj = []
    cr = []
    # add the mutation component
    for i, (AB, Ab, aB, ab) in enumerate(M):
        if AB > 0:
            _coo_helper(ci, cj, cr, i, T[AB-1, Ab+1, aB,   ab  ], alpha*AB)
            _coo_helper(ci, cj, cr, i, T[AB-1, Ab,   aB+1, ab  ], alpha*AB)
        if Ab > 0:
            _coo_helper(ci, cj, cr, i, T[AB+1, Ab-1, aB,   ab  ], alpha*Ab)
            _coo_helper(ci, cj, cr, i, T[AB,   Ab-1, aB,   ab+1], alpha*Ab)
        if aB > 0:
            _coo_helper(ci, cj, cr, i, T[AB+1, Ab,   aB-1, ab  ], alpha*aB)
            _coo_helper(ci, cj, cr, i, T[AB,   Ab,   aB-1, ab+1], alpha*aB)
        if ab > 0:
            _coo_helper(ci, cj, cr, i, T[AB,   Ab+1, aB,   ab-1], alpha*ab)
            _coo_helper(ci, cj, cr, i, T[AB,   Ab,   aB+1, ab-1], alpha*ab)
    # add the drift component
    for i, (X, Y, Z, W) in enumerate(M):
        if X > 0:
            _coo_helper(ci, cj, cr, i, T[X-1, Y+1, Z,   W  ], X*Y)
            _coo_helper(ci, cj, cr, i, T[X-1, Y,   Z+1, W  ], X*Z)
            _coo_helper(ci, cj, cr, i, T[X-1, Y,   Z,   W+1], X*W)
        if Y > 0:
            _coo_helper(ci, cj, cr, i, T[X+1, Y-1, Z,   W  ], Y*X)
            _coo_helper(ci, cj, cr, i, T[X,   Y-1, Z+1, W  ], Y*Z)
            _coo_helper(ci, cj, cr, i, T[X,   Y-1, Z,   W+1], Y*W)
        if Z > 0:
            _coo_helper(ci, cj, cr, i, T[X+1, Y,   Z-1, W  ], Z*X)
            _coo_helper(ci, cj, cr, i, T[X,   Y+1, Z-1, W  ], Z*Y)
            _coo_helper(ci, cj, cr, i, T[X,   Y,   Z-1, W+1], Z*W)
        if W > 0:
            _coo_helper(ci, cj, cr, i, T[X+1, Y,   Z,   W-1], W*X)
            _coo_helper(ci, cj, cr, i, T[X,   Y+1, Z,   W-1], W*Y)
            _coo_helper(ci, cj, cr, i, T[X,   Y,   Z+1, W-1], W*Z)
    # return the coo_matrix
    return scipy.sparse.coo_matrix(
            (cr, (ci, cj)), (nstates, nstates), dtype=float)


##############################################################################
# Get the response.

def guess_denominator(v):
    """
    @param v: a finite distribution
    @return: the smallest positive integer n so that v*n has integer values
    """
    for n in itertools.count(1):
        if n > 1e6:
            raise Exception('failed to guess a denominator')
        vn = v*n
        if np.allclose(vn, np.round(vn), rtol=1e-10):
            return n

def get_response_content(fs):
    # define some parameters
    alpha = 1.0
    #N = fs.N
    k = 4
    # set up print options
    np.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()
    for N in (1, 2, 3, 4):
        print >> out, 'N:', N
        # do everything else
        M = np.array(list(gen_states(N, k)), dtype=int)
        nstates = M.shape[0]
        T = get_inverse_dict(M)
        R_coo = create_coo_moran(M, T, alpha)
        RT_csc = scipy.sparse.csc_matrix(R_coo.T)
        W, V = scipy.sparse.linalg.eigs(RT_csc, k=1, which='SM')
        v = abs(V[:, 0])
        v /= np.sum(v)
        # guess the denominator
        d = guess_denominator(v)
        numerators = np.array(np.round(d*v), dtype=int)
        rat_err = np.max(np.abs(v - numerators/float(d)))
        # print the stuff
        print >> out, 'max abs rationalization error:', rat_err
        print >> out, 'denominator:', d
        headers = (
                'AB', 'Ab', 'aB', 'ab',
                'numerator',
                )
        print >> out, '\t'.join(headers)
        for i in range(nstates):
            p = numerators[i]
            X, Y, Z, W = M[i]
            row = (X, Y, Z, W, p)
            print >> out, '\t'.join(str(x) for x in row)
        print >> out
    return out.getvalue()

