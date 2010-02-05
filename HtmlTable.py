"""
This is a simple module for making an html table.
"""

import unittest
import cgi
from StringIO import StringIO

import Util

def get_labeled_table_string(column_labels, row_labels, row_major_matrix):
    """
    The upper left corner of the matrix will be empty.
    The elements of the matrix will be turned into cgi escaped strings.
    @param column_labels: a label for each column of the matrix
    @param row_labels: a label for each row of the matrix
    @param row_major_matrix: a row major matrix of things to go in the table
    """
    labeled_matrix = []
    labeled_matrix.append([''] + list(column_labels))
    for row_label, row in zip(row_labels, row_major_matrix):
        labeled_matrix.append([row_label] + list(row))
    return get_table_string(labeled_matrix)

def get_table_string(row_major_matrix):
    """
    The elements of the matrix will be turned into cgi escaped strings.
    This is an extremely basic function.
    @param row_major_matrix: a row major matrix of things to go in the table
    """
    out = StringIO()
    print >> out, '<table border>'
    for row in row_major_matrix:
        print >> out, '<tr>'
        for element in row:
            print >> out, '<td>', cgi.escape(str(element)), '</td>'
        print >> out, '</tr>'
    print >> out, '</table>'
    return out.getvalue().strip()


class TestHtmlTable(unittest.TestCase):

    def test_placeholder(self):
        pass


def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestHtmlTable)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

