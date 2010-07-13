"""Create a matrix representation of a Sierpinski-like graph.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import MatrixUtil

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = [
            Form.Integer('iterations', 'number of iterations',
                2, low=0, high=5),
            Form.RadioGroup('format', 'output options', [
                Form.RadioItem('adjacency', 'adjacency matrix'),
                Form.RadioItem('laplacian', 'laplacian matrix', True)])]
    return form_objects

def get_form_out():
    return FormOut.Matrix()

def add_sierpinski(M, offset, iterations):
    if iterations > 0:
        blocksize = 3**(iterations - 1)
        # add the first link
        xoffset = blocksize
        yoffset = blocksize / 2
        M[offset + xoffset][offset + yoffset] = 1
        # add the middle link
        xoffset = 2*blocksize
        yoffset = blocksize - 1
        M[offset + xoffset][offset + yoffset] = 1
        # add the last link
        xoffset = 2*blocksize + blocksize / 2
        yoffset = 2*blocksize - 1
        M[offset + xoffset][offset + yoffset] = 1
        # add the three sub blocks
        for i in range(3):
            add_sierpinski(M, offset + blocksize*i, iterations-1)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # begin the response
    out = StringIO()
    # create the adjacency matrix
    n = 3**fs.iterations
    M = np.zeros((n, n))
    add_sierpinski(M, 0, fs.iterations)
    M = M + M.T
    if fs.adjacency:
        print >> out, MatrixUtil.m_to_string(M)
    elif fs.laplacian:
        for i, row in enumerate(M):
            M[i][i] = -sum(row)
        M = -M
        print >> out, MatrixUtil.m_to_string(M)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
