"""Compute principal components given a .hud file.

The output is an R frame.
The principal components are computed according to Patterson et al.
in their work on population structure.
This is really only for diploid biallelic data.
Microsatellite data is condoned.
Using this for multi-allelic high-ploidy data is probably a hack.
"""

from StringIO import StringIO
import math
import os
import argparse

import numpy as np

from SnippetUtil import HandlingError
import Carbone
import hud
import EigUtil
import Form
import FormOut
import eigenpop

g_tags = ['pca:compute']


g_default_hud_string = """
IC31 1 1 0 0
IC32 1 1 1 0
IC33 1 0 1 1
IC34 0 0 1 0
""".strip()

"""
Form.RadioGroup('scaling', 'PC scaling', [
    Form.RadioItem('scale_none',
        'use orthonormal eigenvectors', True),
    Form.RadioItem('scale_sqrt',
        'scale by the square root of the eigenvalue'),
    Form.RadioItem('scale_eigenvalue',
        'scale by the eigenvalue')]),
"""

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'contents of a .hud file',
                g_default_hud_string),
            Form.CheckGroup('input_options', 'input options', [
                Form.CheckItem('diploid_and_biallelic',
                    'the data source is really diploid and biallelic', True)]),
            Form.Integer('npcs',
                'find this many principal components', 3),
            Form.CheckGroup('output_options', 'output options', [
                Form.CheckItem('add_indices',
                    'add row indices for R table compatibility', True),
                Form.CheckItem('clean_isolates',
                    'force first-column elements to be IC-prefixed', True)])]
    return form_objects

def get_form_out():
    return FormOut.RTable('out')

def get_response_content(fs):
    return process(fs, fs.hud.splitlines()) + '\n'

def process(args, raw_hud_lines):
    """
    @param args: user options from the web or cmdline
    @param hud_lines: raw lines of a .hud file
    @return: results in convenient text form
    """
    out = StringIO()
    names, data = hud.decode(raw_hud_lines)
    # normalize the names of the isolates
    if args.clean_isolates:
        names = [Carbone.clean_isolate_element(x) for x in names]
    # get the pcs
    C_full = np.array(data, dtype=float)
    pcs = eigenpop.get_scaled_eigenvectors(C_full, args.diploid_and_biallelic)
    # check for sufficient number of eigenvectors
    if len(pcs) < args.npcs:
        msg_a = 'the number of requested principal components '
        msg_b = 'must be no more than the number of OTUs'
        raise ValueError(msg_a + msg_b)
    # create the R frame
    headers = ['otu'] + ['pc%d' % (i+1) for i in range(args.npcs)]
    print >> out, '\t'.join(headers)
    for i, name in enumerate(names):
        typed_row = [name] + [pcs[j][i] for j in range(args.npcs)]
        if args.add_indices:
            typed_row = [i+1] + typed_row
        row = [str(x) for x in typed_row]
        print >> out, '\t'.join(row)
    return out.getvalue()
