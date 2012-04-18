"""Convert a .hud file to a .phy file by squishing the columns together.
"""

from SnippetUtil import HandlingError
import Form
import FormOut
import hud
import Phylip

g_tags = ['pca:convert']

g_lines = [
        'IC31 2 0 0 2 0 0 2 0 0 0 1 1 0 2 0 0 2 0',
        'IC32 0 2 0 0 2 0 0 2 0 0 1 1 0 2 0 0 2 0',
        'IC33 0 0 2 0 0 2 0 0 1 1 0 0 2 0 1 1 0 2']

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'table', '\n'.join(g_lines))]
    return form_objects

def get_form_out():
    return FormOut.Hud('out')

def get_response_content(fs):
    headers, data_rows = hud.decode(fs.table.splitlines())
    sequences = [''.join(str(x) for x in row) for row in data_rows]
    return Phylip.encode(headers, sequences)
