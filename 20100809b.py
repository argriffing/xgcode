"""Compute the allele sharing distance given a diploid .hud table.

The allele sharing distance between two individuals
is averaged over all alleles, and for each allele
it is computed as a distance of 0 if there are two alleles in common,
a distance of 1 if there is one allele in common, and a distance of 2 if
there are no alleles in common.
"""

from StringIO import StringIO
import itertools

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import hud

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

def validate_diploid_data_rows(data_rows):
    observed_elements = set(itertools.chain.from_iterable(data_rows))
    expected_elements = set((0,1,2))
    bad_elements = observed_elements - expected_elements
    if bad_elements:
        msg_a = 'expected diploid counts but found '
        msg_b = '{' + ', '.join(str(x) for x in bad_elements) + '}'
        raise ValueError(msg_a + msg_b)

def get_response_content(fs):
    headers, data_rows = hud.decode(fs.table.splitlines())
    validate_diploid_data_rows(data_rows)
    nheaders = len(headers)
    D = np.zeros((nheaders, nheaders))
    for i in range(nheaders):
        for j in range(nheaders):
            ri = np.array(data_rows[i])
            rj = np.array(data_rows[j])
            D[i,j] = np.mean(np.abs(rj - ri))
    return '\n'.join('\t'.join(str(x) for x in r) for r in D)
