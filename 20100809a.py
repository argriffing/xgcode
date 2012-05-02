"""
Create the transpose of a .hud table.
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

g_a_lines = [
        'IC31 2 0 0 2 0 0 2 0 0 0 1 1 0 2 0 0 2 0',
        'IC32 0 2 0 0 2 0 0 2 0 0 1 1 0 2 0 0 2 0',
        'IC33 0 0 2 0 0 2 0 0 1 1 0 0 2 0 1 1 0 2']

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'table', '\n'.join(g_a_lines))]
    return form_objects

def get_form_out():
    return FormOut.Report('out')

def get_response_content(fs):
    headers, data_rows = hud.decode(fs.table.splitlines())
    data_transpose = zip(*data_rows)
    out = StringIO()
    print >> out, ' '.join(headers)
    for row in data_transpose:
        print >> out, ' '.join(str(x) for x in row)
    return out.getvalue()
