"""Sort the names in a non-interleaved phylip alignment.
"""

from SnippetUtil import HandlingError
import Form
import FormOut
import Phylip
import const

g_tags = ['pca:misc']

g_default_data = const.read('20100625a')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                g_default_data),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Phylip()

def get_response_content(fs):
    headers, sequences = Phylip.decode(fs.phylip.splitlines())
    new_pairs = sorted(zip(headers, sequences))
    new_headers, new_sequences = zip(*new_pairs)
    return Phylip.encode(new_headers, new_sequences) + '\n'
