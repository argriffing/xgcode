"""
Read and write .hud files for the fungus project.
"""

from collections import defaultdict

import Util


class HudError(Exception): pass

def decode(raw_lines):
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
            raise HudError(
                    'all rows of a .hud table '
                    'should have the same number of elements')
    headers = [row[0] for row in rows]
    header_counts = defaultdict(int)
    for h in headers:
        header_counts[h] += 1
    repeats = [k for k, v in header_counts.items() if v > 1]
    if len(repeats) > 5:
        raise HudError('%d repeated OTUs within a table' % len(repeats))
    elif repeats:
        raise HudError('repeated OTUs within a table: ' + ', '.join(repeats))
    data = [row[1:] for row in rows]
    for row in data:
        for element in row:
            if element not in list('012'):
                raise HudError(
                        'invalid diploid or haploid element: ' + element)
    data = [[int(x) for x in row] for row in data]
    return headers, data

def encode(headers, rows):
    """
    @param headers: strings
    @param rows: a row major binary matrix
    @return: the contents of a .hud file except for trailing newline
    """
    if not headers:
        return ''
    if len(headers) != len(rows):
        raise ValueError(
                'expected the number of headers (%d) to be the same '
                'as the number of rows (%d)' % (len(headers), len(rows)))
    n = max(len(x) for x in headers)
    ljust_headers = (x.ljust(n+1) for x in headers)
    data_lines = [' '.join(str(x) for x in row) for row in rows]
    return '\n'.join(h + str(x) for h, x in zip(ljust_headers, data_lines))
