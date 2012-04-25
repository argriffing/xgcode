"""
These functions are used directly by the web script framework.
"""

from StringIO import StringIO
import unittest

class MatrixIOError(Exception): pass

def m_to_string(matrix, precision=None):
    """
    @param matrix: a numpy array or row major list of lists
    @param precision: None or the precision for the format conversion
    @return: a multi-line string
    """
    if precision is None:
        return '\n'.join('\t'.join(str(element) for element in row) for row in matrix)
    else:
        format_string = '%%.%df' % precision
        return '\n'.join('\t'.join((format_string % v) for v in row) for row in matrix)

def read_matrix(lines):
    """
    This can be used to extract a matrix from a cgi field or a text file.
    @param lines: lines representing rows of the matrix
    @return: a row major list of lists of floating point elements
    """
    matrix = []
    if not lines:
        raise MatrixIOError('no matrix was specified')
    for i, line in enumerate(lines):
        line = line.strip()
        if line:
            row = []
            for element_string in line.split():
                try:
                    element = float(element_string)
                except ValueError as e:
                    error_lines = [
                            'error on line ' + str(i),
                            line,
                            'this element could not be interpreted as a number: ' + str(element_string)
                            ]
                    raise MatrixIOError('\n'.join(error_lines))
                row.append(element)
            matrix.append(row)
    if not matrix:
        raise MatrixIOError('the matrix is empty')
    row_length_set = set(len(row) for row in matrix)
    if len(row_length_set) != 1:
        raise MatrixIOError('all rows should have the same number of elements')
    return matrix


class TestMatrixIO(unittest.TestCase):

    def test_read_matrix_good(self):
        s = '1 2 3 \n 4 5 6'
        row_major_observed = read_matrix(StringIO(s))
        row_major_expected = [[1, 2, 3], [4, 5, 6]]
        self.assertEquals(row_major_observed, row_major_expected)

    def test_read_matrix_bad(self):
        sio = StringIO('1 2 3 \n 4 5')
        self.assertRaises(MatrixIOError, read_matrix, sio)


if __name__ == '__main__':
    unittest.main()
