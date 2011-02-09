"""Check boundaries of eigenvectors of a simple tree with studded edges.

If the integer branch lengths are positive,
the tree graph has degree sequence 3,2,....,2,1,1,1.
All edge distances and weights will be 1,
at least in the first version.
To estimate first and second derivatives of eigenvectors on a branch,
the branch must have length at least one or two respectively.
"""

import math
from StringIO import StringIO

import numpy as np

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
            Form.Integer('eigk', 'show this eigenvector (Fiedler is 1)',
                1, low=1),
            #Form.Integer('branchk', 'show this branch (first is 1)',
            #   1, low=1, high=3),
            Form.CheckGroup('mycheck', 'extra output options', [
                Form.CheckItem('showmatrix', 'show Laplacian', False)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

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
    # define the total distance of the constructed tree
    d = float(N-1)
    h = 1/d
    # check input compatibility
    if not (fs.eigk+1 <= N):
        msg_a = 'attempting to find a too highly indexed eigenvector '
        msg_b = 'for the number of vertices in the graph'
        raise ValueError(msg_a + msg_b)
    # construct the studded tree Laplacian matrix
    L = create_laplacian_matrix(fs.lena, fs.lenb, fs.lenc)
    # compute the eigendecomposition
    ws, vs = EigUtil.eigh(L)
    # reorder the eigenvalues and eigenvectors
    ws = ws[:-1][::-1]
    vs = vs[:-1][::-1]
    scaling_factor = math.sqrt(N * 0.5)
    # get the eigenvector of interest
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
    if fs.showmatrix:
        print >> out, 'Laplacian matrix:'
        print >> out, L
        print >> out
    return out.getvalue()
