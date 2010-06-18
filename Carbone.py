"""
Add stuff specific to this PCA fungus project.
"""

import re

import numpy as np

import Util

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
        if len(column) != len(set(column)):
            msg = 'expected the column to have unique elements'
            raise RTableError(msg)
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
