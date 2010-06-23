"""
Add stuff specific to this PCA fungus project.
"""

import re
from collections import defaultdict
import itertools

import numpy as np

import Util
import iterutils


g_header_pattern = r'^[a-zA-Z][a-zA-Z0-9]*$'

def is_valid_header(h):
    return re.match(g_header_pattern, h)


def clean_isolate_element(ic):
    if ic is None:
        return ic
    elif ic.startswith('IC'):
        return ic
    else:
        return 'IC' + ic

def clean_isolate_row(row):
    """
    Clean the first element of the row.
    @param row: list whose elements are strings or None
    @return: a new row
    """
    if not row:
        return []
    return [clean_isolate_element(row[0])] + row[1:]

def clean_isolate_table(table):
    """
    Force elements of the first column to start with 'IC'.
    @param table: list of lists whose elements are strings or None
    @return: a new table
    """
    return [clean_isolate_row(row) for row in table]


class RTableError(Exception): pass

class RTable(object):
    def __init__(self, raw_lines):
        """
        Initialize some member variables.
        The member .headers is a list of header strings.
        The member .data is a row major list of list of string data elements
        The member .h_to_i maps a header string to a column index
        """
        self._parse_r_table(raw_lines)
        self.h_to_i = dict((h, i+1) for i, h in enumerate(self.headers))

    def _parse_r_table(self, raw_lines):
        """
        Parse and R table into a header row and data rows.
        Each element of the data table is a string.
        Error checking is minimal.
        """
        lines = Util.get_stripped_lines(raw_lines)
        header_line, data_lines = lines[0], lines[1:]
        self.headers = header_line.split()
        self.data = [line.split() for line in data_lines]
        nheaders = len(self.headers)
        if len(set(self.headers)) != nheaders:
            msg = 'multiple columns are labeled with the same header'
            raise RTableError(msg)
        for row in self.data:
            if len(row) != nheaders+1:
                msg_a = 'the header row has %d elements ' % nheaders
                msg_b = 'and a data row has %d elements; ' % len(row)
                msg_c = 'all data rows should have one more element '
                msg_d = 'than the header row'
                raise RTableError(msg_a + msg_b + msg_c + msg_d)
    
    def header_to_column_index(self, header):
        if header not in self.h_to_i:
            msg = 'the column header %s was not found in the table' % header
            raise RTableError(msg)
        return self.h_to_i[header]

    def header_to_column(self, header):
        column_index = self.header_to_column_index(header)
        column = [row[column_index] for row in self.data]
        return column

    def header_to_primary_column(self, header):
        column_index = self.header_to_column_index(header)
        column = [row[column_index] for row in self.data]
        d = defaultdict(int)
        for x in column:
            d[x] += 1
        repeated_keys = [k for k, v in d.items() if v > 1]
        if len(repeated_keys) > 5:
            msg_a = '%d repeated keys ' % len(repeated_keys)
            msg_b = 'found in the primary column'
            raise RTableError(msg_a + msg_b)
        elif repeated_keys:
            msg_a = 'repeated keys in the primary column: '
            msg_b = ', '.join(repeated_keys)
            raise RTableError(msg_a + msg_b)
        return column


class WordError(Exception): pass

class Word(object):
    def __init__(self, line):
        # store the original line
        self.line = line
        # extract the name and the binary vector
        elements = line.split()
        self.name = elements[0]
        for x in elements[1:]:
            if x not in ('0', '1'):
                msg = 'expected 0 or 1 but got ' + x
                raise WordError(msg)
        self.v = np.array([int(x) for x in elements[1:]])

def validate_words(words):
    """
    Make sure the binary vectors are the same lengths.
    @param words: word objects
    """
    lengths = set(len(w.v) for w in words)
    if len(lengths) != 1:
        msg = 'each binary vector should be the same length'
        raise WordError(msg)

def get_words(lines):
    """
    @param lines: lines of a .hud file
    @return: a list of validated words
    """
    lines = Util.get_stripped_lines(lines)
    words = [Word(line) for line in lines]
    validate_words(words)
    return words

class HudError(Exception): pass

def read_hud(raw_lines):
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


def get_compound_column(column, uniques):
    """
    This is a step in converting a multivalued matrix to a binary matrix.
    Note that 'uniques' may have elements not in 'column'.
    @param column: a list
    @param uniques: a list of ordered unique elements at this position
    @return: a list of binary lists
    """
    n = len(uniques)
    element_to_index = dict((x, i) for i, x in enumerate(uniques))
    index_column = [element_to_index[x] for x in column]
    compound_column = []
    for index in index_column:
        sublist = [(1 if i == index else 0) for i in range(n)]
        compound_column.append(sublist)
    return compound_column

def get_binary_rows_helper(columns, uniques):
    """
    Convert multivalued data to binary data.
    @param multivalued_columns: list of lists whose elements are anything
    @param uniques: list of unique element lists
    @return: longer rows of binary elements
    """
    compound_columns = [
            get_compound_column(x, u) for x, u in zip(columns, uniques)]
    # convert to compound binary rows
    compound_rows = zip(*compound_columns)
    # convert to simple binary rows
    return [list(itertools.chain(*x)) for x in compound_rows]

def get_binary_rows(multivalued_rows):
    """
    Convert multivalued data to binary data.
    @param multivalued_rows: elements of each rows can be anything
    @return: longer rows of binary elements
    """
    # get the columns
    columns = zip(*multivalued_rows)
    # convert to compound binary columns
    uniques = [list(iterutils.unique_everseen(x)) for x in columns]
    # convert to simple binary rows
    return get_binary_rows_helper(columns, uniques)

def get_hud_content(headers, rows):
    """
    @param headers: strings
    @param rows: a row major binary matrix
    @return: the content of a .hud file
    """
    n = max(len(x) for x in headers)
    ljust_headers = (x.ljust(n+1) for x in headers)
    data_lines = [' '.join(str(x) for x in row) for row in rows]
    return '\n'.join(h + str(x) for h, x in zip(ljust_headers, data_lines))
