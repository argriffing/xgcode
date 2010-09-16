"""Rank loci by their correlation with a given principal component. [unfinished]
"""

from StringIO import StringIO
import math
import os

import numpy as np
import argparse

from SnippetUtil import HandlingError
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

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'contents of a .hud file', g_default_hud_string),
            Form.CheckGroup('input_options', 'input options', [
                Form.CheckItem('diploid_and_biallelic',
                    'the data source is really diploid and biallelic', True)]),
            Form.Integer('axis',
                'rank by correlation with this axis (the first axis is 1)', 1),
            Form.RadioGroup('locus_options', 'locus indexing output format', [
                Form.RadioItem('locus_from_0', 'count loci from 0', True),
                Form.RadioItem('locus_from_1', 'count loci from 1')])
            Form.RadioGroup('rank_options', 'ranking output format', [
                Form.RadioItem('rank_squared',
                    'rank by decreasing squared correlation', True),
                Form.RadioItem('rank_unsquared',
                    'rank by decreasing correlation')])
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    return process(fs, fs.hud.splitlines()) + '\n'

def mycorr(a, b):
    stda = np.std(a)
    stdb = np.std(b)
    if not stda:
        return 0
    if not stdb:
        return 0
    m = np.mean((a - np.mean(a)) * (b - np.mean(b)))
    return m / (stda * stdb)

def process(args, raw_hud_lines):
    """
    @param args: user options from the web or cmdline
    @param hud_lines: raw lines of a .hud file
    @return: results in convenient text form
    """
    out = StringIO()
    names, data = hud.decode(raw_hud_lines)
    pcs = eigenpop.get_scaled_eigenvectors(data, args.diploid_and_biallelic)
    # check for sufficient number of eigenvectors
    if len(evecs) < args.ncoords:
        msg_a = 'the number of requested principal components '
        msg_b = 'must be no more than the number of OTUs'
        raise ValueError(msg_a + msg_b)
    # compute the correlation of each SNP vector with each principal PC
    mylist = []
    for snp in C_full.T:
        row = [mycorr(snp, pc) for pc in pcs[:args.ncoords]]
        mylist.append(row)
    np.set_printoptions(linewidth=300, threshold=10000)
    return str(np.array(mylist))
