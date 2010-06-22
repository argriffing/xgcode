"""Convert a phylip alignment to a binary character alignment.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import Phylip
import Carbone


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

def process(raw_lines):
    headers, sequences = Phylip.read_non_interleaved(raw_lines)
    binary_rows = Carbone.get_binary_rows(sequences)
    return Carbone.get_hud_content(headers, binary_rows)
