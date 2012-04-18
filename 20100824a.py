"""Update a combination log for a .phy to .hud conversion.

The combination log file looks really free-form so I
am probably making unrealistic assumptions about its format.
I'll assume that the paragraph after the 'Loci combined' line
is the relevant one, and I will ignore the taxon list.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import Phylip
import Carbone
import Util
import iterutils
import hud
import const

g_tags = ['pca:misc']

g_phy_lines = [
        '4 4',
        'IC31      AC1C',
        'IC32      AC2C',
        'IC33      AT3G',
        'IC34      AT4T']

g_combo_lines = [
        'Loci combined',
        '',
        'aflM_GlobalSpecies_alignment.phy  1-1',
        'aflW_GlobalSpecies_alignment.phy  2-2',
        'amdS_GlobalSpecies_alignment.phy  3-4',
        '',
        '',
        'Taxon list',
        '',
        '1 IC31',
        '2 IC32',
        '3 IC33',
        '',
        '4 IC34']


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('phylip', 'non-interleaved Phylip alignment',
                '\n'.join(g_phy_lines)),
            Form.MultiLine('combo', 'locus combination log file',
                '\n'.join(g_combo_lines))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def gen_combo_line_triples(raw_lines):
    """
    Yield (name, low, high) triples where low and high are 1-based.
    @param raw_lines: lines of a locus combination file
    """
    paragraphs = list(Util.gen_paragraphs(raw_lines))
    if paragraphs[0] != ['Loci combined']:
        msg_a = "expected the only line of the first paragraph "
        msg_b = "of the combination input to be 'Loci combined'"
        raise ValueError(msg_a + msg_b)
    for line in paragraphs[1]:
        elements = line.split()
        if len(elements) != 2:
            msg_a = 'expected two whitespace separated elements '
            msg_b = 'on this line of the combination input: '
            msg_c = line
            raise ValueError(msg_a + msg_b + msg_c)
        name, lowhigh = elements
        if lowhigh.count('-') != 1:
            msg_a = 'expected a single hyphen in the second '
            msg_b = 'whitespace separated element '
            msg_c = 'on this line of the combination input: '
            msg_d = line
            raise ValueError(msg_a + msg_b + msg_c + msg_d)
        s_low, s_high = lowhigh.split('-')
        try:
            low = int(s_low)
            high = int(s_high)
        except ValueError as e:
            msg = 'expected integer ranges in the combination input'
            raise ValueError(msg)
        if low < 1 or high < 1:
            msg_a = 'expected positive integer range bounds '
            msg_b = 'in the combination input'
            raise ValueError(msg_a + msg_b)
        if high < low:
            msg_a = 'the lower bound of a combination range '
            msg_b = 'should not be greater than the upper bound'
            raise ValueError(msg_a + msg_b)
        yield name, low, high

def expand_ranges(ranges):
    """
    Blindly yield elements inside the inclusive ranges.
    @param low_high_pairs: a sequence of (low, high) inclusive ranges
    """
    for low, high in low_high_pairs:
        for j in range(low, high+1):
            yield j

def get_response_content(fs):
    # get the combo info
    combo_triples = list(gen_combo_line_triples(fs.combo.splitlines()))
    names, lows, highs = zip(*combo_triples)
    ranges = zip(lows, highs)
    if lows[0] != 1:
        msg = 'expected the first lower bound to be 1'
        raise ValueError(msg)
    for (low, high), (nlow, nhigh) in iterutils.pairwise(ranges):
        if high + 1 != nlow:
            msg_a = 'expected the next lower bound '
            msg_b = 'to be one more than the current upper bound'
            raise ValueError(msg_a + msg_b)
    # get the phylip info
    headers, sequences = Phylip.decode(fs.phylip.splitlines())
    phylip_columns = zip(*sequences)
    counts = [len(set(col)) for col in phylip_columns]
    # validate the compatibility between the combo and phylip data
    if highs[-1] != len(phylip_columns):
        msg_a = 'expected the last upper bound to be '
        msg_b = 'equal to the number of columns of the phylip alignment'
        raise ValueError(msg_a + msg_b)
    # get the sum of counts in each combination group
    combo_counts = []
    for i, (low, high) in enumerate(ranges):
        combo_count = 0
        # note that low and high are 1-based and inclusive
        for j in range(low-1, high):
            combo_count += counts[j]
        combo_counts.append(combo_count)
    # write the new combo log
    out = StringIO()
    print >> out, 'Loci combined'
    print >> out
    k = 0
    for name, count in zip(names, combo_counts):
        low = k + 1
        high = k + count
        print >> out, '%s\t%d-%d' % (name, low, high)
        k += count
    return out.getvalue()
