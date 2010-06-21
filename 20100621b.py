"""Convert a phylip alignment to a binary character alignment.
"""

from StringIO import StringIO
import itertools

from SnippetUtil import HandlingError
import Form
import Util
import Phylip
import iterutils


g_default_lines = [
        '3 20',
        'Cow       ATGGCATATCCCATACAACT',
        'Carp      ATGGCACACCCAACGCAACT',
        'Chicken   ATGGCCAACCACTCCCAACT']

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                '\n'.join(g_default_lines)),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.phylip.splitlines())
    filename = 'out.hud'
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def get_compound_column(column):
    """
    @param column: a list
    @return: a list of binary lists
    """
    uniq = list(iterutils.unique_everseen(column))
    n = len(uniq)
    element_to_index = dict((x, i) for i, x in enumerate(uniq))
    index_column = [element_to_index[x] for x in column]
    compound_column = []
    for index in index_column:
        sublist = [(1 if i == index else 0) for i in range(n)]
        compound_column.append(sublist)
    return compound_column

def process(raw_lines):
    # read the headers and the sequences
    headers, sequences = Phylip.read_non_interleaved(raw_lines)
    # get the columns
    columns = zip(*sequences)
    # convert to compound binary columns
    compound_columns = [get_compound_column(x) for x in columns]
    # convert to compound binary rows
    compound_rows = zip(*compound_columns)
    # convert to simple binary rows
    binary_rows = [list(itertools.chain(*x)) for x in compound_rows]
    # create the .hud content
    n = max(len(x) for x in headers)
    ljust_headers = (x.ljust(n+1) for x in headers)
    data_lines = [' '.join(str(x) for x in row) for row in binary_rows]
    return '\n'.join(h + str(x) for h, x in zip(ljust_headers, data_lines))
