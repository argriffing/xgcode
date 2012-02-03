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

def validate_headers(headers):
    for h in headers:
        if not is_valid_header(h):
            raise ValueError('invalid column header: %s' % h)

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


class NumericError(Exception): pass

def get_numeric_column(data, index):
    """
    @param data: row major list of lists of numbers as strings
    @param index: column index
    @return: list of floats
    """
    strings = zip(*data)[index]
    try:
        floats = [float(x) for x in strings]
    except ValueError, v:
        raise NumericError
    return floats

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
