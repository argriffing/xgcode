"""
Look at properties of Hessians computed by algopy for molecular models.

One question is how various orders of Pade-based expm approximations compare to
the eigh-based expm approximations for small time-reversible models,
with respect to accuracy of gradients and Hessians.
Perhaps this depends on the distinctness of the eigenvalues
or the order of the Pade approximation.
In this script, avoid the numpy abbreviation np to help make the code
look more like the algopy code.
The molecular model in this toy example is HKY85.
"""

from StringIO import StringIO

import numpy
import scipy
import algopy

import Form
import FormOut

# AGCT subsititution counts between human and chimp mitochondrial coding dna.
g_data = numpy.array([
        [2954, 141, 17, 16],
        [165, 1110, 5, 2],
        [18, 4, 3163, 374],
        [15, 2, 310, 2411],
        ], dtype=float)

def algopy_expm(A, q):
    """
    Compute the matrix exponential using Pade approximation.

    Reference
    ---------

    N. J. Higham,
    "The Scaling and Squaring Method for the Matrix Exponential Revisited",
    SIAM. J. Matrix Anal. & Appl. 26, 1179 (2005).
    """
    if not (A.ndim == 2 and A.shape[0] == A.shape[1]):
        raise ValueError('expected a square matrix')
    N = A.shape[0]
    ident = numpy.eye(N)
    if q == 7:
        pade = _algopy_pade7
    elif q == 13:
        pade = _algopy_pade13
    else:
        raise NotImplementedError('see recent publications by Nick Higham')
    U, V = pade(A, ident)
    return algopy.solve(-U + V, U + V)

def _algopy_pade7(A, ident):
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    A2 = algopy.dot(A, A)
    A4 = algopy.dot(A2, A2)
    A6 = algopy.dot(A2, A4)
    U = algopy.dot(A, b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    return U, V

def _algopy_pade13(A, ident):
    b = (
            64764752532480000., 32382376266240000., 7771770303897600.,
            1187353796428800., 129060195264000., 10559470521600.,
            670442572800., 33522128640., 1323241920.,
            40840800., 960960., 16380., 182., 1.)
    A2 = algopy.dot(A, A)
    A4 = algopy.dot(A2, A2)
    A6 = algopy.dot(A2, A4)
    U = algopy.dot(A,
            algopy.dot(A6, b[13]*A6 + b[11]*A4 + b[9]*A2) + (
                b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident))
    V = algopy.dot(A6, b[12]*A6 + b[10]*A4 + b[8]*A2) + (
            b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident)
    return U, V

def transform_params(Y):
    X = algopy.exp(Y)
    tsrate, tvrate = X[0], X[1]
    v_unnormalized = algopy.zeros(4, dtype=X)
    v_unnormalized[0] = X[2]
    v_unnormalized[1] = X[3]
    v_unnormalized[2] = X[4]
    v_unnormalized[3] = 1.0
    v = v_unnormalized / algopy.sum(v_unnormalized)
    return tsrate, tvrate, v

def get_numpy_rate_matrix(Y):
    a, b, v = transform_params(Y)
    Q = numpy.array([
        [0, a, b, b],
        [a, 0, b, b],
        [b, b, 0, a],
        [b, b, a, 0],
        ])
    Q = Q * v
    Q -= numpy.diag(numpy.sum(Q, axis=1))
    return Q

def eval_f_expm_scipy(Y):
    a, b, v = transform_params(Y)
    Q = numpy.array([
        [0, a, b, b],
        [a, 0, b, b],
        [b, b, 0, a],
        [b, b, a, 0],
        ])
    Q = Q * v
    Q -= numpy.diag(numpy.sum(Q, axis=1))
    P = scipy.linalg.expm(Q)
    S = numpy.log(numpy.dot(numpy.diag(v), P))
    return -numpy.sum(S * g_data)

def eval_f_expm(Y, pade_order):
    a, b, v = transform_params(Y)
    Q = algopy.zeros((4,4), dtype=Y)
    Q[0,0] = 0;    Q[0,1] = a;    Q[0,2] = b;    Q[0,3] = b;
    Q[1,0] = a;    Q[1,1] = 0;    Q[1,2] = b;    Q[1,3] = b;
    Q[2,0] = b;    Q[2,1] = b;    Q[2,2] = 0;    Q[2,3] = a;
    Q[3,0] = b;    Q[3,1] = b;    Q[3,2] = a;    Q[3,3] = 0;
    Q = Q * v
    Q -= algopy.diag(algopy.sum(Q, axis=1))
    P = algopy_expm(Q, pade_order)
    S = algopy.log(algopy.dot(algopy.diag(v), P))
    return -algopy.sum(S * g_data)

def eval_f_eigh(Y):
    a, b, v = transform_params(Y)
    Q = algopy.zeros((4,4), dtype=Y)
    Q[0,0] = 0;    Q[0,1] = a;    Q[0,2] = b;    Q[0,3] = b;
    Q[1,0] = a;    Q[1,1] = 0;    Q[1,2] = b;    Q[1,3] = b;
    Q[2,0] = b;    Q[2,1] = b;    Q[2,2] = 0;    Q[2,3] = a;
    Q[3,0] = b;    Q[3,1] = b;    Q[3,2] = a;    Q[3,3] = 0;
    Q = algopy.dot(Q, algopy.diag(v))
    Q -= algopy.diag(algopy.sum(Q, axis=1))
    va = algopy.diag(algopy.sqrt(v))
    vb = algopy.diag(1./algopy.sqrt(v))
    W, U = algopy.eigh(algopy.dot(algopy.dot(va, Q), vb))
    M = algopy.dot(U, algopy.dot(algopy.diag(algopy.exp(W)), U.T))
    P = algopy.dot(vb, algopy.dot(M, va))
    S = algopy.log(algopy.dot(algopy.diag(v), P))
    return -algopy.sum(S * g_data)

def eval_grad(Y, fn, *moreargs):
    """
    compute the gradient of fn in the forward mode of AD
    """
    Y = algopy.UTPM.init_jacobian(Y)
    retval = fn(Y, *moreargs)
    return algopy.UTPM.extract_jacobian(retval)

def eval_hess(Y, fn, *moreargs):
    """
    compute the hessian of fn in the forward mode of AD
    """
    Y = algopy.UTPM.init_hessian(Y)
    retval = fn(Y, *moreargs)
    return algopy.UTPM.extract_hessian(Y.shape[0], retval)

def eval_grad_eigh(Y):
    return eval_grad(Y, eval_f_eigh)

def eval_hess_eigh(Y):
    return eval_hess(Y, eval_f_eigh)

def eval_grad_expm(Y, q):
    return eval_grad(Y, eval_f_expm, q)

def eval_hess_expm(Y, q):
    return eval_hess(Y, eval_f_expm, q)



def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.RadioGroup('pade_options', 'expm Pade approximation', [
                Form.RadioItem('q7', 'order 7', True),
                Form.RadioItem('q13', 'order 13')]),
            Form.RadioGroup('point_options', 'get Taylor info at this point', [
                Form.RadioItem('use_guess', 'initial guess', True),
                Form.RadioItem('use_mle', 'maximum likelihood estimate')]),
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    numpy.set_printoptions(linewidth=200)
    out = StringIO()
    if fs.q7:
        q = 7
    elif fs.q13:
        q = 13
    else:
        raise ValueError
    Y = numpy.zeros(5)
    #
    if fs.use_mle:
        results = scipy.optimize.fmin_ncg(
            eval_f_eigh,
            Y,
            eval_grad_eigh,
            fhess_p=None,
            fhess=eval_hess_eigh,
            args=(),
            avextol=1e-05,
            epsilon=1.4901161193847656e-08,
            maxiter=100,
            full_output=True,
            disp=1,
            retall=0,
            callback=None)
        Y = results[0]
    #
    Q = get_numpy_rate_matrix(Y)
    W = scipy.linalg.eigvals(Q)
    #
    print >> out, '--------------------------------'
    print >> out, 'properties of the function and the point'
    print >> out, 'at which its Taylor expansion is evaluated'
    print >> out, '--------------------------------'
    print >> out, 'Q\n', Q
    print >> out, 'Q - Q.T\n', Q - Q.T
    print >> out, 'eigenvalues of Q\n', W
    print >> out, 'Y\n', Y
    print >> out, '--------------------------------'
    print >> out
    print >> out, '--------------------------------'
    print >> out, ' simple check (functions)       '
    print >> out, '--------------------------------'
    print >> out, 'eval_f_expm_scipy(Y)\n', eval_f_expm_scipy(Y)
    print >> out, 'eval_f_expm(Y, q)\n', eval_f_expm(Y, q)
    print >> out, 'eval_f_eigh(Y)\n', eval_f_eigh(Y)
    print >> out, 'eval_f_eigh(Y) - eval_f_expm(Y, q)\n', (
            eval_f_eigh(Y) - eval_f_expm(Y, q))
    print >> out, 'eval_f_eigh(Y) - eval_f_expm_scipy(Y)\n', (
            eval_f_eigh(Y) - eval_f_expm_scipy(Y))
    print >> out, 'eval_f_expm(Y, q) - eval_f_expm_scipy(Y)\n', (
            eval_f_expm(Y, q) - eval_f_expm_scipy(Y))
    print >> out, '--------------------------------'
    print >> out
    print >> out, '--------------------------------'
    print >> out, ' simple check (gradients)       '
    print >> out, '--------------------------------'
    print >> out, 'eval_grad_expm(Y, q)\n', eval_grad_expm(Y, q)
    print >> out, 'eval_grad_eigh(Y)\n', eval_grad_eigh(Y)
    print >> out, 'eval_grad_eigh(Y) - eval_grad_expm(Y, q)\n', (
            eval_grad_eigh(Y) - eval_grad_expm(Y, q))
    print >> out, '--------------------------------'
    print >> out
    print >> out, '--------------------------------'
    print >> out, ' simple check (hessians)        '
    print >> out, '--------------------------------'
    print >> out, 'eval_hess_expm(Y, q)\n', eval_hess_expm(Y, q)
    print >> out, 'eval_hess_eigh(Y)\n', eval_hess_eigh(Y)
    print >> out, 'eval_hess_eigh(Y) - eval_hess_expm(Y, q)\n', (
            eval_hess_eigh(Y) - eval_hess_expm(Y, q))
    print >> out, '--------------------------------'
    print >> out
    #
    return out.getvalue()

