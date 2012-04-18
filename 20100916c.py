"""Rank loci by their correlation with a given principal component.
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
            Form.Integer('axis',
                'rank by correlation with this axis (the first axis is 1)', 1),
            Form.CheckGroup('input_options', 'input options', [
                Form.CheckItem('diploid_and_biallelic',
                    'the data source is really diploid and biallelic', True)]),
            Form.RadioGroup('locus_options', 'locus indexing output format', [
                Form.RadioItem('locus_from_0', 'count loci from 0', True),
                Form.RadioItem('locus_from_1', 'count loci from 1')]),
            Form.RadioGroup('rank_options', 'ranking output format', [
                Form.RadioItem('rank_squared',
                    'rank by decreasing squared correlation', True),
                Form.RadioItem('rank_unsquared',
                    'rank by decreasing correlation')])]
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
    C_full = np.array(data, dtype=float)
    pcs = eigenpop.get_scaled_eigenvectors(C_full, args.diploid_and_biallelic)
    axis_index = args.axis - 1
    # check for sufficient number of eigenvectors
    if axis_index >= len(pcs):
        msg = 'the requested axis is not available'
        raise ValueError(msg)
    # compute the correlation of each SNP vector the requested PC
    pc = pcs[axis_index]
    corrs = [mycorr(snp, pc) for snp in C_full.T]
    sqcorrs = [mycorr(snp, pc)**2 for snp in C_full.T]
    if args.rank_squared:
        keys = sqcorrs
    else:
        keys = corrs
    corr_index_pairs = [(cor, i) for i, cor in enumerate(keys)]
    sorted_pairs = list(reversed(sorted(corr_index_pairs)))
    indices = zip(*sorted_pairs)[1]
    if args.locus_from_1:
        nominal_indices = [i+1 for i in indices]
    else:
        nominal_indices = indices
    rows = [(nom_i, corrs[i]) for i, nom_i in zip(indices, nominal_indices)]
    lines = ['\t'.join(str(x) for x in row) for row in rows]
    return '\n'.join(lines) + '\n'
