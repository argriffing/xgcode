"""Combine binary character sequence copies.

Sequences whose names are the same except
for a terminal lowercase letter will be combined using the selected operation.
The input .hud table should reflect binary characters
or possibly counts of binary characters in the case of non-haploid data.
"""

import collections

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import hud
import iterutils

g_tags = ['carbone_lab']

g_lines = [
        'IC31 2 0 0 2 0 0 2 0 0 0 1 1 0 2 0 0 2 0',
        'IC32 0 2 0 0 2 0 0 2 0 0 1 1 0 2 0 0 2 0',
        'IC33 0 0 2 0 0 2 0 0 1 1 0 0 2 0 1 1 0 2']


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud', 'binary character table',
                '\n'.join(g_lines)),
            Form.RadioGroup('options', 'sequence combination operation', [
                Form.RadioItem('combine_exist',
                    'combine by existence', True),
                Form.RadioItem('combine_count',
                    'combine by counting')]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report('out.hud', [])

def process_headers(headers):
    """
    Detect header variants.
    Try to stably group the headers.
    Each element in the output list is a pair.
    The first element of the pair is the header prefix.
    The second element of the pair is a set of whole headers.
    @param headers: sequence of headers
    @return: a reordered, reduced, and augmented list
    """
    prefixes = []
    for h in headers:
        prefix = None
        if len(h) > 1 and h[-1].islower():
            prefix = h[:-1]
        prefixes.append(prefix)
    # detect prefixes that occur more than once
    precount = collections.defaultdict(int)
    for p in prefixes:
        if p:
            precount[p] += 1
    common_prefixes = set(k for k, v in precount.items() if v > 1)
    # map frequent prefixes to header sets
    prefix_to_headers = collections.defaultdict(set)
    for h, p in zip(headers, prefixes):
        if p in common_prefixes:
            prefix_to_headers[p].add(h)
    # get the final data structure
    used_prefixes = set()
    pairs = []
    for h, p in zip(headers, prefixes):
        if p in used_prefixes:
            continue
        if p in common_prefixes:
            pairs.append([p, prefix_to_headers[p]])
            used_prefixes.add(p)
        else:
            pairs.append([h, [h]])
    return pairs

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the headers and data from all of the input sources
    headers, sequences = hud.decode(fs.hud.splitlines())
    h_to_s = dict((h, s) for h, s in zip(headers, sequences))
    headers_out = []
    sequences_out = []
    for p, hs in process_headers(headers):
        headers_out.append(p)
        data = np.vstack(h_to_s[h] for h in hs).sum(axis=0)
        if fs.combine_exist:
            data = np.minimum(1, data)
        sequences_out.append(data)
    text = hud.encode(headers_out, sequences_out) + '\n'
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'out.hud') 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, text
