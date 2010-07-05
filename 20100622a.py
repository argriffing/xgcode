"""Convert diploid microsatellite data to a ternary character alignment.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Form
import Carbone
import hud
import Util
import iterutils

g_default_rows = [
        ['DOE3-0002Aa', '94', '94', '156', '162', '89', '172'],
        ['DOE3-0002Ab', '94', '94', '156', '182', '89', '172'],
        ['DOE3-0002Ba', '91', '91', '?', '162', '89', '172'],
        ['DOE3-0002Bb', '91', '91', '?', '182', '89', '172'],
        ['DOE3-0007Aa', '89', '89', '158', '160', '93', '?'],
        ['DOE3-0007Ab', '89', '89', '160', '160', '107', '?']]

g_default_lines = ['\t'.join(x for x in row) for row in g_default_rows]

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('data', 'diploid microsatellite data',
                '\n'.join(g_default_lines)),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.data.splitlines())
    filename = 'out.hud'
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def read_satellite_lines(raw_lines):
    """
    How can i combine the two haploid data sources?
    Maybe create each data matrix separately from the interleaved input.
    @param raw_lines: raw input lines
    @return: headers, diploid data
    """
    lines = Util.get_stripped_lines(raw_lines)
    if len(lines) % 2:
        raise ValueError('expected an even number of lines')
    if len(lines) < 2:
        raise ValueError('expected at least two lines')
    full_rows = [x.split() for x in lines]
    nfullcols = len(full_rows[0])
    if nfullcols < 2:
        raise ValueError('expected at least two columns')
    for row in full_rows:
        if len(row) != nfullcols:
            msg = 'each row should have the same number of elements'
            raise ValueError(msg)
    a_full_rows = [row for i, row in enumerate(full_rows) if i % 2 == 0]
    b_full_rows = [row for i, row in enumerate(full_rows) if i % 2 == 1]
    a_headers = [row[0] for row in a_full_rows]
    b_headers = [row[0] for row in b_full_rows]
    for h in a_headers:
        if not h.endswith('a'):
            msg = 'each odd row label should end with the letter a'
            raise ValueError(msg)
    for h in b_headers:
        if not h.endswith('b'):
            msg = 'each even row label should end with the letter b'
            raise ValueError(msg)
    headers = [h[:-1] for h in a_headers]
    # get the unique elements of each column
    rows = [row[1:] for row in full_rows]
    cols = zip(*rows)
    uniques = [list(iterutils.unique_everseen(col)) for col in cols]
    # get the results for each row
    a_rows = [row[1:] for row in a_full_rows]
    b_rows = [row[1:] for row in b_full_rows]
    a_columns = zip(*a_rows)
    b_columns = zip(*b_rows)
    a_binary_rows = Carbone.get_binary_rows_helper(a_columns, uniques)
    b_binary_rows = Carbone.get_binary_rows_helper(b_columns, uniques)
    # add the elements entrywise and return as a list of lists
    binary_rows = (np.array(a_binary_rows) + np.array(b_binary_rows)).tolist()
    return headers, binary_rows

def process(raw_lines):
    headers, binary_rows = read_satellite_lines(raw_lines)
    return hud.encode(headers, binary_rows) + '\n'
