"""
Utility functions for interfacing with R.
"""

import unittest

def get_table_string(M, column_headers):
    """
    Convert a row major rate matrix to a string representing an R table.
    @param M: a row major matrix
    @param column_headers: the labels of the data columns
    """
    # assert that the matrix is rectangular
    assert len(set(len(row) for row in M)) == 1
    # assert that the number of columns is the same as the number of column headers
    assert len(M[0]) == len(column_headers)
    # assert that the headers are valid
    for header in column_headers:
        if '_' in header:
            raise ValueError('the header "%s" is invalid because it uses an underscore' % header)
    # define each line in the output string
    lines = []
    lines.append('\t'.join([''] + list(column_headers)))
    for i, row in enumerate(M):
        R_row = [float_to_R(value) for value in row]
        lines.append('\t'.join([str(i+1)] + R_row))
    return '\n'.join(lines)

def float_to_R(value):
    """
    Convert a python floating point value to a string usable by R.
    @param value: a floating point number
    @return: the R string representing the floating point number
    """
    if value == float('inf'):
        return 'Inf'
    else:
        return str(value)

class TestRUtil(unittest.TestCase):

    def test_placeholder(self):
        """
        do nothing
        """
        pass


if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestRUtil)
    unittest.TextTestRunner(verbosity=2).run(suite)


