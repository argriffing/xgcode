"""
Read and write .hud files for the fungus project.
"""

from collections import defaultdict

import Util


class HudError(Exception): pass

def parse(raw_lines):
    """
    @param raw_lines: raw lines of a hud file
    @return: headers, data
    """
    lines = Util.get_stripped_lines(raw_lines)
    if not lines:
        return [], []
    rows = [line.split() for line in lines]
    ncols = len(rows[0])
    for row in rows:
        if len(row) != ncols:
            msg_a = 'all rows of a .hud table '
            msg_b = 'should have the same number of elements'
            raise HudError(msg_a + msg_b)
    headers = [row[0] for row in rows]
    header_counts = defaultdict(int)
    for h in headers:
        header_counts[h] += 1
    repeats = [k for k, v in header_counts.items() if v > 1]
    if len(repeats) > 5:
        msg = '%d repeated OTUs within a table' % len(repeats)
        raise HudError(msg)
    elif repeats:
        msg = 'repeated OTUs within a table: ' + ', '.join(repeats)
        raise HudError(msg)
    data = [row[1:] for row in rows]
    for row in data:
        for element in row:
            if element not in list('012'):
                msg = 'invalid diploid or haploid element: ' + element
                raise HudError(msg)
    data = [[int(x) for x in row] for row in data]
    return headers, data

def to_blob(headers, rows):
    """
    @param headers: strings
    @param rows: a row major binary matrix
    @return: the content of a .hud file
    """
    if not headers:
        return ''
    n = max(len(x) for x in headers)
    ljust_headers = (x.ljust(n+1) for x in headers)
    data_lines = [' '.join(str(x) for x in row) for row in rows]
    return '\n'.join(h + str(x) for h, x in zip(ljust_headers, data_lines))
