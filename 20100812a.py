"""Convert diploid SNP data from AWclust format to .hud format.

The value -1 is used for missing data.
"""

from StringIO import StringIO
import os
import itertools

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import hud
import iterutils

g_tags = ['pca:convert']

g_data = """IC31 IC32 IC33
2 0 0
0 2 0
0 0 2
2 0 0
0 2 0
0 0 2
2 0 0
0 2 0
0 -1 -1
0 -1 -1
1 1 0
1 1 0
0 0 2
2 2 0
0 0 1
0 0 1
2 2 0
0 0 2"""


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
        Form.MultiLine('data', 'SNP data in AWclust format', g_data),
        Form.CheckGroup('options', 'missing data options', [
            Form.CheckItem('use_mode',
                'replace missing data with the non-missing mode',
                True)])]
    return form_objects

def get_form_out():
    return FormOut.Report('out')

def get_response_content(fs):
    lines = Util.get_stripped_lines(fs.data.splitlines())
    if len(lines) < 2:
        raise ValueError('expected at least two lines')
    rows = [line.split() for line in lines]
    headers = rows[0]
    data_rows = [[int(x) for x in row] for row in rows[1:]]
    for row in data_rows:
        for x in row:
            if x not in (-1, 0, 1, 2):
                msg = 'invalid diploid data value: %d' % x
                raise ValueError(msg)
    # impute the missing data
    if fs.use_mode:
        imputed_data_rows = []
        for row in data_rows:
            non_missing_row = [x for x in row if x != -1]
            if not non_missing_row:
                msg = 'a variable has missing data for each individual'
                raise ValueError(msg)
            counts = [0]*3
            for x in non_missing_row:
                counts[x] += 1
            imputed_value = counts.index(max(counts))
            imputed_row = [imputed_value if x == -1 else x for x in row]
            imputed_data_rows.append(imputed_row)
        data_rows = imputed_data_rows
    # return the hud table
    return hud.encode(headers, zip(*data_rows))
