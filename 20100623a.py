"""Combine .hud tables.

A .hud table is an alignment of binary (haploid) or ternary (diploid) traits.
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

g_tags = ['ztools:convert']

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

def get_form_out():
    return FormOut.Report('out')

def get_response_content(fs):
    return process([fs.table_a.splitlines(), fs.table_b.splitlines()]) + '\n'

def process(line_sources):
    """
    @param line_sources: sources of line iterables
    """
    # get the headers and data from all of the input sources
    header_data_pairs = [hud.decode(lines) for lines in line_sources]
    header_list, data_list = zip(*header_data_pairs)
    # get the header to index map for each input source
    h_to_i_list = [Util.inverse_map(x) for x in header_list]
    # get the intersection of headers in all lists
    header_sets = [set(x) for x in header_list]
    header_intersection = set.intersection(*header_sets)
    # get the ordered list of all headers
    unique_headers = list(iterutils.unique_everseen(
            itertools.chain.from_iterable(header_list)))
    # get the ordered list of headers present in every input source
    out_headers = [h for h in unique_headers if h in header_intersection]
    out_data = []
    for h in out_headers:
        row = []
        for data, h_to_i in zip(data_list, h_to_i_list):
            if h in h_to_i:
                row.extend(data[h_to_i[h]])
        out_data.append(row)
    return hud.encode(out_headers, out_data) + '\n'

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infiles', type=argparse.FileType('r'), nargs='+',
            help='.hud files to be combined')
    args = parser.parse_args()
    main(args)
    for f in args.infiles:
        f.close()
