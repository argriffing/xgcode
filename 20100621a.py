"""Add row labels to a table to satisfy the R table requirements.

Input should be a table of whitespace separated elements with column headers.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import Util

g_tags = ['carbone_lab']


g_default_rows = [
        ['temperature', 'precipitation', 'otu'],
        ['15', '600', 'IC31'],
        ['30', '700', 'IC32'],
        ['45', '800', 'IC33']]

g_default_lines = ['\t'.join(row) for row in g_default_rows]


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R-like table without row labels',
                '\n'.join(g_default_lines)),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable('out.table', [])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs, fs.table.splitlines())
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'out.table') 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def process(args, raw_lines):
    lines = Util.get_stripped_lines(raw_lines)
    rows = [x.split() for x in lines]
    if len(rows) < 2:
        raise ValueError('expected at least two rows')
    ncols = len(rows[0])
    for row in rows:
        if len(row) != ncols:
            raise ValueError('each row should have the same number of columns')
    headers, data = rows[0], rows[1:]
    # define the output header and data
    out_headers = headers
    out_data = [[i+1]+row for i, row in enumerate(data)]
    # write the R table
    out = StringIO()
    print >> out, '\t'.join(out_headers)
    for row in out_data:
        print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue()
