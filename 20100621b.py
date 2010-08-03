"""Convert a phylip alignment to a binary character alignment.
"""

from SnippetUtil import HandlingError
import Form
import FormOut
import Phylip
import Carbone
import hud
import const

g_tags = ['carbone_lab']

g_default_data = const.read('20100625a')

#FIXME multiple output types

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                g_default_data),
            Form.RadioGroup('format_out', 'output format', [
                Form.RadioItem('hud', 'hud', True),
                Form.RadioItem('phy', 'non-interleaved phylip')]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report('out.%s', None, ['format_out'])

def get_response_content(fs):
    return process(fs, fs.phylip.splitlines()) + '\n'

def process(fs, raw_lines):
    headers, sequences = Phylip.decode(raw_lines)
    binary_rows = Carbone.get_binary_rows(sequences)
    if fs.hud:
        return hud.encode(headers, binary_rows) + '\n'
    elif fs.phy:
        binary_seqs = [''.join(str(x) for x in row) for row in binary_rows]
        return Phylip.encode(headers, binary_seqs) + '\n'
