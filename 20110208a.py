"""
Check boundaries of eigenvectors of a simple tree with studded edges.

If the integer branch lengths are positive,
the tree graph has degree sequence 3,2,....,2,1,1,1.
All edge distances and weights will be 1,
at least in the first version.
To estimate first and second derivatives of eigenvectors on a branch,
the branch must have length at least one or two respectively.
Apparently the scipy sparse eigendecomposition is buggy.
So if you use a sparse matrix representation
then take the result with a grain of salt.
Edited to not use an obsolete scipy function name
http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=621056
"""

import math
from StringIO import StringIO

import numpy as np
import scipy.sparse
import scipy.sparse.linalg

import Form
import FormOut
import Euclid
import EigUtil
import iterutils

def get_form():
    """
    @return: a list of form objects
    """
    # define the form objects
    form_objects = [
            Form.Integer('lena', 'integer length of first branch',
                5, low=0, high=400),
            Form.Integer('lenb', 'integer length of second branch',
                7, low=0, high=400),
            Form.Integer('lenc', 'integer length of third branch',
                13, low=0, high=400),
            Form.Integer('eigk', 'use this eigenvector (Fiedler is 1)',
                1, low=1),
            #Form.Integer('branchk', 'show this branch (first is 1)',
            #   1, low=1, high=3),
            Form.CheckGroup('mycheck', 'extra options', [
                Form.CheckItem('sparse', 'use a sparse matrix format', True),
                Form.CheckItem('showv', 'show the eigenpair', False),
                Form.CheckItem('showmatrix', 'show the matrix', False)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def create_laplacian_csr_matrix(lena, lenb, lenc):
    """
    Create a sparse matrix.
    Begin by creating the triangular matrix where nonzero elements
    have row index less than column index.
    Then add the opposite triangle to form a symmetric matrix.
    Then add the diagonal elements.
    @param lena: integer length of first branch.
    @param lenb: integer length of second branch.
    @param lenc: integer length of third branch.
    """
    # Initialize the triples: matrix element, row index, and column index.
    drc_triples = []
    # Initialize the data list, row index list, and column index list.
    #d, r, c = []
    # Initialize the number of vertices.
    N = 1 + lena + lenb + lenc
    # Add connections to the hub vertex.
    if lena:
        drc_triples.append((-1, 0, 1))
    if lenb:
        drc_triples.append((-1, 0, lena+1))
    if lenc:
        drc_triples.append((-1, 0, lena+lenb+1))
    # Add connections on the first branch.
    for i in range(lena-1):
        j = i + 1
        drc_triples.append((-1, j, j+1))
    # Add connections on the second branch.
    for i in range(lenb-1):
        j = lena + i + 1
        drc_triples.append((-1, j, j+1))
    # Add connections on the second branch.
    for i in range(lenc-1):
        j = lena + lenb + i + 1
        drc_triples.append((-1, j, j+1))
    # Define the opposite triangle drc triples.
    drc_opp = [(d, c, r) for d, r, c in drc_triples]
    # Define the diagonal drc triples.
    drc_diag_a = [(-d, c, c) for d, r, c in drc_triples]
    drc_diag_b = [(-d, c, c) for d, r, c in drc_opp]
    # Define the total drc triples.
    drc_total = drc_triples + drc_opp + drc_diag_a + drc_diag_b
    # Define the separate {d, r, c} lists.
    vd, vr, vc = zip(*drc_total)
    # Create the sparse coordinate format matrix.
    S = scipy.sparse.coo_matrix((np.array(vd, dtype=float), (vr, vc)), (N, N))
    # Return a more sophisticated compressed sparse row format matrix.
    return S.tocsr()

def create_laplacian_matrix(lena, lenb, lenc):
    """
    @param lena: integer length of first branch.
    @param lenb: integer length of second branch.
    @param lenc: integer length of third branch.
    """
    N = 1 + lena + lenb + lenc
    A = np.zeros((N,N), dtype=float)
    # Add connections to the hub vertex.
    if lena:
        A[0, 1] = 1
        A[1, 0] = 1
    if lenb:
        A[0, lena+1] = 1
        A[lena+1, 0] = 1
    if lenc:
        A[0, lena+lenb+1] = 1
        A[lena+lenb+1, 0] = 1
    # Add tridiagonal connections on the first branch.
    for i in range(lena-1):
        j = i + 1
        A[j, j+1] = 1
        A[j+1, j] = 1
    # Add tridiagonal connections on the second branch.
    for i in range(lenb-1):
        j = lena + i + 1
        A[j, j+1] = 1
        A[j+1, j] = 1
    # Add tridiagonal connections on the second branch.
    for i in range(lenc-1):
        j = lena + lenb + i + 1
        A[j, j+1] = 1
        A[j+1, j] = 1
    L = Euclid.adjacency_to_laplacian(A)
    return L

class BranchInfo:
    def __init__(self):
        self.p0 = None
        self.p1 = None
        self.p2 = None
        self.q0 = None
        self.q1 = None
        self.q2 = None
        self.width = None
        self.k = None

def value_to_string(x):
    if x is None:
        return 'undefined'
    return str(x)

def get_response_content(fs):
    # define the number of nodes
    N = 1 + fs.lena + fs.lenb + fs.lenc
    # check input compatibility
    if not (fs.eigk+1 <= N):
        msg_a = 'attempting to find a too highly indexed eigenvector '
        msg_b = 'for the number of vertices in the graph'
        raise ValueError(msg_a + msg_b)
    if N < 2:
        raise ValueError('the tree has no length')
    # define the total distance of the constructed tree
    d = float(N-1)
    h = 1/d
    # construct the studded tree Laplacian matrix
    if fs.sparse:
        v0 = np.ones(N, dtype=float)
        L_csr = create_laplacian_csr_matrix(fs.lena, fs.lenb, fs.lenc)
        arpack_k = fs.eigk+1
        ncv = 3*arpack_k + 3
        ws, vs = scipy.sparse.linalg.eigsh(
                L_csr, arpack_k, which='SM',
                v0=v0,
                ncv=ncv, return_eigenvectors=True)
        ws = ws[1:]
        vs = vs.T[1:]
    else:
        L = create_laplacian_matrix(fs.lena, fs.lenb, fs.lenc)
        ws, vs = EigUtil.eigh(L)
        ws = ws[:-1][::-1]
        vs = vs[:-1][::-1]
    scaling_factor = math.sqrt(N * 0.5)
    # get the eigenvector of interest
    eigenvalue = ws[fs.eigk-1]
    v = vs[fs.eigk-1]
    # init the branch info
    binfos = [BranchInfo() for i in range(3)]
    for i, binfo in enumerate(binfos):
        binfo.k = i+1
        # split the eigenvector of interest into the branch components
        if binfo.k == 1:
            offset = 1
            binfo.width = fs.lena
            w = np.array([v[0]] + v[offset:offset+binfo.width].tolist())
        elif binfo.k == 2:
            offset = 1 + fs.lena
            binfo.width = fs.lenb
            w = np.array([v[0]] + v[offset:offset+binfo.width].tolist())
        elif binfo.k == 3:
            offset = 1 + fs.lena + fs.lenb
            binfo.width = fs.lenc
            w = np.array([v[0]] + v[offset:offset+binfo.width].tolist())
        else:
            raise ValueError
        # compute some boundary info
        if len(w) >= 1:
            binfo.p0 = w[0]
        if len(w) >= 2:
            binfo.p1 = (w[1] - w[0]) / h
        if len(w) >= 3:
            binfo.p2 = (w[0] - 2*w[1] + w[2]) / (h*h)
        if len(w) >= 1:
            binfo.q0 = w[-1]
        if len(w) >= 2:
            binfo.q1 = (w[-1] - w[-2]) / h
        if len(w) >= 3:
            binfo.q2 = (w[-3] - 2*w[-2] + w[-1]) / (h*h)
    # begin writing the report
    np.set_printoptions(linewidth=200, threshold=10000)
    out = StringIO()
    # summarize global properties
    print >> out, 'total branch length:'
    print >> out, N - 1
    print >> out
    print >> out, 'total number of graph vertices including degree 2 vertices:'
    print >> out, N
    print >> out
    # show the sum of first derivatives near the hub
    if N > 1:
        p1sum = 0
        for binfo in binfos:
            if binfo.p1:
                p1sum += binfo.p1
        p1sum_string = str(p1sum)
    else:
        d1sum_string = 'undefined'
    print >> out, "sum of f'(x) on all branches near the hub:", p1sum_string
    print >> out
    # summarize properties per branch per eigenvector
    for binfo in binfos:
        print >> out, 'summary of eigenvector', fs.eigk, 'on branch', binfo.k
        print >> out, 'unscaled branch length:', binfo.width
        if binfo.width:
            print >> out, 'internal', ''.join(['-']*binfo.width), 'pendant'
            print >> out, "internal f(x):  ", value_to_string(binfo.p0)
            print >> out, "internal f'(x): ", value_to_string(binfo.p1)
            print >> out, "internal f''(x):", value_to_string(binfo.p2)
            print >> out, "pendant  f(x):  ", value_to_string(binfo.q0)
            print >> out, "pendant  f'(x): ", value_to_string(binfo.q1)
            print >> out, "pendant  f''(x):", value_to_string(binfo.q2)
        print >> out
    if fs.showv:
        print >> out, 'the eigenvalue:'
        print >> out, eigenvalue
        print >> out
        print >> out, 'the whole eigenvector:'
        print >> out, v
        print >> out
    if fs.showmatrix:
        if fs.sparse:
            print >> out, 'Laplacian matrix (from sparse internal repr):'
            print >> out, L_csr.toarray()
            print >> out
        else:
            print >> out, 'Laplacian matrix (from dense internal repr):'
            print >> out, L
            print >> out
    return out.getvalue()
