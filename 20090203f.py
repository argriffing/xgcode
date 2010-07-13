"""Given a Laplacian matrix, use principal minors to find the distance matrix.

Given a Laplacian matrix,
use ratios of principal minors to find the distance matrix.
This uses an incorrect formula from "On Euclidean distance matrices"
and a corrected version by Eric Stone.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    L = np.array([
            [1.0, 0.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, -1.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, -1.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, -1.0],
            [-1.0, -1.0, 0.0, 0.0, 3.0, -1.0],
            [0.0, 0.0, -1.0, -1.0, -1.0, 3.0]])
    # define the form objects
    form_objects = [
            Form.Matrix('laplacian', 'combinatorial Laplacian matrix',
                L, MatrixUtil.assert_symmetric),
            Form.RadioGroup('method', 'method choices', [
                Form.RadioItem('incorrect', 'use the method from the paper'),
                Form.RadioItem('correct', 'use the corrected method', True)])]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def get_deleted_matrix(M, row_indices, column_indices):
    """
    @param M: a numpy array
    @param row_indices: the indices of the rows to delete
    @param column_indices: the indices of the columns to delete
    """
    nrows, ncols = M.shape
    D = []
    for i in range(nrows):
        row = []
        for j in range(ncols):
            if j not in column_indices:
                row.append(M[i][j])
        if i not in row_indices:
            D.append(row)
    return np.array(D)

def get_minor(M, row_indices, column_indices):
    """
    Get the determinant of the matrix with rows and columns deleted.
    @param M: a square matrix
    @param row_indices: the indices of the rows to delete
    @param column_indices: the indices of the columns to delete
    """
    M_deleted = get_deleted_matrix(M, row_indices, column_indices)
    return np.linalg.det(M_deleted)

def get_incorrect_distance_matrix(L):
    """
    Use the formula from "On Euclidean distance matrices".
    @param L: a Laplacian matrix
    """
    n = len(L)
    D = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i][j] = get_minor(L, [i], [j]) / get_minor(L, [i], [i])
    return D

def get_correct_distance_matrix(L):
    """
    Use a corrected formula.
    @param L: a Laplacian matrix
    """
    n = len(L)
    D = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i][j] = get_minor(L, [i, j], [i, j]) / get_minor(L, [i], [i])
    return D

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix
    L = fs.laplacian
    # get the distance matrix
    if fs.correct:
        D = get_correct_distance_matrix(L)
    else:
        D = get_incorrect_distance_matrix(L)
    # define the response
    out = StringIO()
    print >> out, MatrixUtil.m_to_string(D)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
