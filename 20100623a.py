"""Combine .hud tables.

A .hud table is an alignment of binary (haploid) or ternary (diploid) traits.
"""

from StringIO import StringIO
import os

from SnippetUtil import HandlingError
import Form
import Util
import Carbone
import iterutils

g_a_lines = [
        'IC31 2 0 0 2 0 0 2 0 0 0 1 1 0 2 0 0 2 0',
        'IC32 0 2 0 0 2 0 0 2 0 0 1 1 0 2 0 0 2 0',
        'IC33 0 0 2 0 0 2 0 0 1 1 0 0 2 0 1 1 0 2']

g_b_lines = [
        'IC31 1 1 0 0',
        'IC32 1 1 1 0',
        'ICXX 1 1 1 1',
        'IC33 1 0 1 1']


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table_a', 'first .hud table',
                '\n'.join(g_a_lines)),
            Form.MultiLine('table_b', 'second .hud table',
                '\n'.join(g_b_lines)),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs, fs.table_a.splitlines(), fs.table_b.splitlines())
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'out.hud') 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text

def process(args, raw_a_lines, raw_b_lines):
    a_headers, a_data = Carbone.read_hud(raw_a_lines)
    b_headers, b_data = Carbone.read_hud(raw_b_lines)
    a_h_to_i = dict((h, i) for i, h in enumerate(a_headers))
    b_h_to_i = dict((h, i) for i, h in enumerate(b_headers))
    b_set = set(b_headers)
    out_headers = [h for h in a_headers if h in b_set]
    out_data = []
    for h in out_headers:
        row = []
        if h in a_h_to_i:
            row.extend(a_data[a_h_to_i[h]])
        if h in b_h_to_i:
            row.extend(b_data[b_h_to_i[h]])
        out_data.append(row)
    return Carbone.get_hud_content(out_headers, out_data)
