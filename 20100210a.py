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
import JC69


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
    p_ref_change = JC69.distance_to_probability(fs.ref_length)
    p_child_change = JC69.distance_to_probability(fs.child_length)
    # For now sum over all possibilities of non-reference nodes.
    # This could be done more efficiently using Felsenstein pruning,
    # but I am ignoring this for now.
    p_RR = 0.0
    p_RA = 0.0
    p_AA = 0.0
    p_AB = 0.0
    ref = 0
    for c12 in range(4):
        if c12 == ref:
            p12 = 1.0 - p_ref_change
        else:
            p12 = p_ref_change / 3.0
        for c1 in range(4):
            if c1 == c12:
                p1 = p12 * (1.0 - p_child_change)
            else:
                p1 = p12 * (p_child_change / 3.0)
            for c2 in range(4):
                if c2 == c12:
                    p2 = p1 * (1.0 - p_child_change)
                else:
                    p2 = p1 * (p_child_change / 3.0)
                # Classify the joint distribution
                # and add weight to the appropriate state.
                if c1 == ref and c2 == ref:
                    p_RR += p2
                elif c1 == ref or c2 == ref:
                    p_RA += p2
                elif c1 == c2:
                    p_AA += p2
                else:
                    p_AB += p2
    total = p_RR + p_RA + p_AA + p_AB
    if abs(total - 1) > 1e-7:
        raise HandlingError('internal error -- probs do not sum to one')
    # write the distribution
    out = StringIO()
    print >> out, 'p(RR):', p_RR
    print >> out, 'p(RA):', p_RA
    print >> out, 'p(AA):', p_AA
    print >> out, 'p(AB):', p_AB
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
