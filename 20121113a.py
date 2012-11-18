"""
Check 2-parameter recessivity estimates with gtr mutation and cond. log lik.
"""

from StringIO import StringIO
import os
import collections

import numpy

import Form
import FormOut
import npcodon
import yangdata

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()

def get_selection_F(log_counts, compo, log_nt_weights):
    return log_counts - numpy.dot(compo, log_nt_weights)

def get_selection_S(F):
    e = numpy.ones_like(F)
    return numpy.outer(e, F) - numpy.outer(F, e)

def get_D(
        gtr, compo,
        log_counts,
        d, log_kb, log_nt_weights, log_repop,
        ):
    #
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F) * numpy.exp(log_repop)
    soft_sign_S = numpy.tanh(numpy.exp(log_kb)*S)
    D = d * soft_sign_S
    return D

def construct_args():
    Args = collections.namedtuple('Args',
            ['fmin', 'disease', 'mtdna', 'infile', 't1', 't2'])
    args = Args(
            fmin = None,
            disease = 'kacser',
            mtdna = None,
            infile = os.path.expanduser(
                '~/data/YangNielsen2008MBE.MutSel/HCMMR.txt'),
            t1 = 'mm8',
            t2 = 'rn4',
            )
    return args

def get_response_content(fs):
    numpy.set_printoptions(
            linewidth=1000000,
            threshold=1000000,
            )
    out = StringIO()
    #
    args = construct_args()
    #
    #
    # Precompute some ndarrays
    # according to properties of DNA and the genetic code.
    if args.mtdna:
        code = npcodon.g_code_mito
        stop = npcodon.g_stop_mito
    else:
        code = npcodon.g_code
        stop = npcodon.g_stop
    #
    all_codons = npcodon.enum_codons(stop)
    codons = all_codons[:-len(stop)]
    gtr = npcodon.get_gtr(codons)
    syn, nonsyn = npcodon.get_syn_nonsyn(code, codons)
    compo = npcodon.get_compo(codons)
    asym_compo = npcodon.get_asym_compo(codons)
    ham = npcodon.get_hamming(codons)
    #
    subs_counts = yangdata.get_subs_counts_from_data_files(args)
    codon_counts = (
            numpy.sum(subs_counts, axis=0) + numpy.sum(subs_counts, axis=1))
    codon_counts = codon_counts[:len(codons)]
    subs_counts = subs_counts[:len(codons), :len(codons)]
    v = codon_counts / float(numpy.sum(codon_counts))
    log_counts = numpy.log(codon_counts)
    #
    log_mu = -4.02576875,
    log_gtr_exch = numpy.array([
        0.50335873,
        1.31231415,
        0.3491126,
        -0.17953527,
        1.55231821,
        0])
    log_omega = -2.2685352
    #
    d = 2.86523675
    log_kb = -0.63087496
    log_nt_weights = numpy.array([
        -0.12604391,
        0.46524165,
        0.96822465,
        0])
    log_repop = -0.01399715
    #
    D = get_D(
            gtr, compo,
            log_counts, d, log_kb, log_nt_weights, log_repop)
    print >> out, D
    print >> out, numpy.unique(D)
    mask = numpy.sum(gtr, axis=2)
    D = D * mask
    print >> out, 'unique in mask:', numpy.unique(mask)
    print >> out, D
    print >> out, numpy.unique(D)
    #
    return out.getvalue()

