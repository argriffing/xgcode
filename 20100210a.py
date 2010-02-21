"""Compute a distribution implied by Jukes Cantor on a simple tree.

The tree is so simple that it has three branches and an internal node.
Two of the branches have the same length.
The tip associated with the third branch is the reference.
Given the two branch lengths that define the tree,
compute the joint distribution over the non-reference tip values.
That is, the distribution over the events {RR, RA, AA, AB}.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import DGRP

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('ref_length', 'length of the reference branch',
                '0.1', low_inclusive=0.0),
            Form.Float('child_length', 'length of each child branch',
                '0.1', low_inclusive=0.0)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    distn = DGRP.get_zygosity_distribution(fs.ref_length, fs.child_length)
    out = StringIO()
    print >> out, 'p(RR):', distn[0]
    print >> out, 'p(RA):', distn[1]
    print >> out, 'p(AA):', distn[2]
    print >> out, 'p(AB):', distn[3]
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
