"""Compute zygosity distributions given four branch length parameters.

The labels of the four parameters are not precise
and should be interpreted only roughly.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import DGRP

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('w', 'effect of misalignment',
                '0.6', low_inclusive=0.0),
            Form.Float('x', 'distance from the reference',
                '0.1', low_inclusive=0.0),
            Form.Float('y', 'distance between lines',
                '0.01', low_inclusive=0.0),
            Form.Float('z', 'distance from the nominal line',
                '0.0001', low_inclusive=0.0)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def distn_to_string(distn):
    names = ('RR', 'RA', 'AA', 'AB')
    lines = ['p(%s): %s' % (name, p) for name, p in zip(names, distn)]
    return '\n'.join(lines)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # heterozygous
    print >> out, 'heterozygous region:'
    a = fs.x
    b = fs.y + fs.z
    distn = DGRP.get_zygosity_distribution(a, b)
    print >> out, distn_to_string(distn)
    print >> out
    # homozygous
    print >> out, 'homozygous region:'
    a = fs.x + fs.y
    b = fs.z
    distn = DGRP.get_zygosity_distribution(a, b)
    print >> out, distn_to_string(distn)
    print >> out
    # misaligned
    print >> out, 'misaligned homozygous region:'
    a = fs.w + fs.x + fs.y
    b = fs.z
    distn = DGRP.get_zygosity_distribution(a, b)
    print >> out, distn_to_string(distn)
    print >> out
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
