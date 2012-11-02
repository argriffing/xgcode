"""
Check selection coefficients estimated within a mutation-selection model.

This currently uses hardcoded max liklihood parameter estimates
from a model with two parameters for recessivity as a function of selection.
The idea is to check the implicit selection estimates as an intermediate step
towards checking the recessivity values implied by the combination of the
recessivity parameter estimates with the estimated selection coefficients.
"""

from StringIO import StringIO
import math
import os
import collections

import numpy
import scipy
import scipy.special

import Form
import FormOut
import kimrecessive
import npcodon
import yangdata

# These are hardcoded max likelihood estimates for an 8 parameter model.
# They control:
# rate, kappa, omega, d, kb, nt_mut_distn, nt_mut_distn, nt_mut_distn
# where kappa is a transition/transversion ratio,
# omega accounts for amino acid selection,
# and d and kb are two parameters that control recessivity.
# Some of the parameters are log transformed.
# These estimates are for the mouse-rat pair in the Nielsen-Yang 2008 data.
g_opt_x = numpy.array(
        [
            -3.69603727,
            1.27847313,
            -2.24987323,
            6.47323139,
            -1.5383088,
            -0.14216812,
            0.37080688,
            0.49925123,
          ], dtype=float)

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

def get_sparse_D(
        ts, tv, compo,
        log_counts,
        d, log_kb, log_nt_weights,
        ):
    #
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F)
    soft_sign_S = numpy.tanh(numpy.exp(log_kb)*S)
    D = d * soft_sign_S
    return (ts + tv) * D

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
    ts, tv = npcodon.get_ts_tv(codons)
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
    log_mu, log_kappa, log_omega, d, log_kb, nta, ntc, ntg = g_opt_x.tolist()
    log_nt_weights = numpy.array([nta, ntc, ntg, 0])
    D = get_sparse_D(
            ts, tv, compo, log_counts, d, log_kb, log_nt_weights)
    print >> out, D
    print >> out, numpy.unique(D)
    #
    return out.getvalue()

