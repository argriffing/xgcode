"""Convert k-ploid microsatellite data to a (k+1)-ary character alignment.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Carbone
import hud
import Util
import iterutils
import const

g_tags = ['pca:misc']

g_default_data = const.read('20100730w')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('data', 'k-ploid microsatellite data',
                g_default_data)]
    return form_objects

def get_form_out():
    return FormOut.Hud('out.hud')

def get_response_content(fs):
    headers, binary_rows = read_microsatellite_lines(fs.data.splitlines())
    return hud.encode(headers, binary_rows) + '\n'

def get_ploidy(names):
    """
    Get the ploidy from the name sequence.
    The sequence of names should be ordered,
    and each name should end with a lower case letter such as 'a'
    or 'b' for diploid data, or 'a', 'b', 'c', 'd' for tetraploid data.
    @param names: the sequence of names.
    """
    ploidy = None
    for i, name in enumerate(names):
        indicator = name[-1]
        offset = ord(indicator) - ord('a')
        if ploidy:
            expected_offset = i % ploidy
            if offset == expected_offset:
                # the indicator is consistent with the ploidy
                continue
            else:
                # the indicator is inconsistent with the ploidy
                msg_a = 'autodetection found a ploidy of %d ' % ploidy
                msg_b = 'so the label of row %d was expected ' % (i+1)
                msg_c = 'to end with "%c"' % (expected_offset + ord('a'))
                raise ValueError(msg_a + msg_b + msg_c)
        elif offset == i:
            # the ploidy is in the process of being detected
            continue
        elif i == 0:
            # the first line is wrong
            raise ValueError('the label of the first line should end with "a"')
        elif offset == 0:
            # the ploidy has been detected
            ploidy = i
            continue
        else:
            # found an error before the ploidy was detected
            lines = [
                    'the label of row %d was expected to end ' % (i+1),
                    'with "a" if the ploidy ',
                    'of the organism is %d ' % i,
                    'or to end with "%c" if the ploidy ' % (ord('a') + i),
                    'of the organism is greater than %d ' % i]
            raise ValueError(' '.join(lines))
    return ploidy

def gen_headers(full_names, ploidy):
    """
    Yield headers without the lowercase suffix.
    @param full_names: a sequence of names including the lowercase suffix
    @param ploidy: a small integer ploidy
    """
    nchunks, remainder = divmod(len(full_names), ploidy)
    if remainder:
        msg = 'the number of rows should be a multiple of the detected ploidy'
        raise ValueError(msg)
    for i in range(nchunks):
        expected_base = full_names[i*ploidy][:-1]
        for j in range(ploidy):
            observed_base = full_names[i*ploidy + j][:-1]
            if observed_base != expected_base:
                msg_a = 'all but the last letter of each row label '
                msg_b = 'should be consistent within each group'
                raise ValueError(msg_a + msg_b)
        yield expected_base

def read_microsatellite_lines(raw_lines):
    """
    How can i combine the two haploid data sources?
    Maybe create each data matrix separately from the interleaved input.
    @param raw_lines: raw input lines
    @return: headers, diploid data
    """
    lines = Util.get_stripped_lines(raw_lines)
    full_rows = [line.split() for line in lines]
    nfullcols = len(full_rows[0])
    if nfullcols < 2:
        raise ValueError('expected at least two columns')
    if not all(len(row) == nfullcols for row in full_rows):
        raise ValueError('expected the same number of elements in each row')
    full_cols = zip(*full_rows)
    full_names = full_cols[0]
    ploidy = get_ploidy(full_names)
    headers = list(gen_headers(full_names, ploidy))
    # get the unique elements of each column
    rows = [row[1:] for row in full_rows]
    cols = zip(*rows)
    uniques = [list(iterutils.unique_everseen(col)) for col in cols]
    # get the rows for each offset
    n = len(rows) / ploidy
    groups = [[rows[j*ploidy + i] for j in range(n)] for i in range(ploidy)]
    # get the column groups
    col_groups = [zip(*m) for m in groups]
    # get the binary row groups
    bin_row_groups = [
            Carbone.get_binary_rows_helper(
                cols, uniques) for cols in col_groups]
    # get the entrywise sum
    binary_rows = np.array(bin_row_groups).sum(axis=0).tolist()
    return headers, binary_rows
