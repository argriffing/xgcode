"""Get a combination log for a .phy to .hud conversion.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import Phylip
import Carbone
import hud
import const

g_phy_lines = [
        '4 4',
        'IC31      AC1C',
        'IC32      AC2C',
        'IC33      AT3G',
        'IC34      AT4T']

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                '\n'.join(g_phy_lines)),
            Form.CheckGroup('out_checkbox', 'uninformative loci', [
                Form.CheckItem('remove', 'remove uninformative loci', True)]),
            Form.RadioGroup('out_radio', 'indexing', [
                Form.RadioItem('from_zero', 'count columns from 0'),
                Form.RadioItem('from_one', 'count columns from 1', True)]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    headers, sequences = Phylip.decode(fs.phylip.splitlines())
    phylip_columns = zip(*sequences)
    counts = [len(set(col)) for col in phylip_columns]
    header_low_high = []
    k = 0 if fs.from_zero else 1
    for i, count in enumerate(counts):
        if not (count == 1 and fs.remove):
            low = k
            high = k + (count - 1)
            header_low_high.append((i+(1 if fs.from_one else 0), low, high))
        k += count
    out = StringIO()
    for triple in header_low_high:
        print >> out, 'locus%s.phy\t%d-%d' % triple
    return out.getvalue()
