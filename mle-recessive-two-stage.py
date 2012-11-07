"""
This is intended to be a faster version of the similarly named script.

It is an exercise in premature optimization.
Or maybe it is not premature.
"""

import string
import math
import argparse
import functools
import warnings

import numpy
import scipy
import scipy.optimize
import scipy.special
import scipy.linalg
import algopy
import algopy.special

import jeffopt
import kimrecessive
import npcodon
import yangdata


##########################################################################
# algopy stuff involving parameters
#
# These two functions could possibly be moved out of this module.
# But they could depend on algopy in other contexts.

def get_selection_F(log_counts, compo, log_nt_weights):
    """
    The F and S notation is from Yang and Nielsen 2008.
    Note that three of the four log nt weights are free parameters.
    One of the four log weights is zero and the other three
    are free parameters to be estimated jointly in the
    maximimum likelihood search,
    so this function is inside the optimization loop.
    @param log_counts: logs of empirical codon counts
    @param compo: codon composition as defined in the get_compo function
    @param log_nt_weights: un-normalized log mutation process probabilities
    @return: a log selection for each codon, up to an additive constant
    """
    return log_counts - numpy.dot(compo, log_nt_weights)

def get_selection_S(F):
    """
    The F and S notation is from Yang and Nielsen 2008.
    @param F: a selection value for each codon, up to an additive constant
    @return: selection differences F_j - F_i, also known as S_ij
    """
    e = numpy.ones_like(F)
    return numpy.outer(e, F) - numpy.outer(F, e)


##########################################################################
# algopy stuff involving parameters

def get_fixation_unconstrained(S, d):
    sign_S = algopy.sign(S)
    D = d * sign_S
    H = algopy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_piecewise(
                    0.5*S[i, j], D[i, j])
    return H

def get_fixation_unconstrained_quad(S, d):
    """
    Use numerical quadrature.
    Do not bother trying to use algopy for this.
    """
    sign_S = numpy.sign(S)
    D = d * sign_S
    H = numpy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_quad(
                    0.5*S[i, j], D[i, j])
    return H

def get_Q_prefix(
        ts, tv, syn, nonsyn,
        log_mu, log_kappa, log_omega):
    """
    Compute a chunk of a hadamard decomposition of the pre-Q matrix.
    By hadamard decomposition I mean the factoring of a matrix
    into the entrywise product of two matrices.
    By pre-Q matrix I mean the rate matrix before the row sums
    have been subtracted from the diagonal.
    Notation is from Yang and Nielsen 2008.
    The first group of args consists of precomputed ndarrays.
    The second group depends only on free parameters.
    Note that this function does not depend on mutation process
    stationary distribution parameters,
    and it does not depend on recessivity parameters.
    """
    mu = algopy.exp(log_mu)
    kappa = algopy.exp(log_kappa)
    omega = algopy.exp(log_omega)
    return mu * (kappa * ts + tv) * (omega * nonsyn + syn)

def get_Q_postfix_recessivity(
        compo, asym_compo,
        h,
        log_counts,
        d, log_nt_weights):
    """
    This is specific to the model of recessivity.
    """
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F)
    pre_Q_postfix = algopy.exp(algopy.dot(asym_compo, log_nt_weights)) * h(S, d)
    return pre_Q_postfix
    Q = pre_Q - algopy.diag(algopy.sum(pre_Q, axis=1))
    return Q

#FIXME: remove
def get_Q_unconstrained(
        ts, tv, syn, nonsyn, compo, asym_compo,
        h,
        log_counts,
        log_mu, log_kappa, log_omega, d, log_nt_weights):
    """
    This adds a single parameter.
    """
    mu = algopy.exp(log_mu)
    kappa = algopy.exp(log_kappa)
    omega = algopy.exp(log_omega)
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F)
    pre_Q = mu * (kappa * ts + tv) * (omega * nonsyn + syn) * algopy.exp(
            algopy.dot(asym_compo, log_nt_weights)) * h(S, d)
    Q = pre_Q - algopy.diag(algopy.sum(pre_Q, axis=1))
    return Q

def get_log_likelihood(P, v, subs_counts):
    """
    The stationary distribution of P is empirically derived.
    It is proportional to the codon counts by construction.
    @param P: a transition matrix using codon counts and free parameters
    @param v: stationary distribution proportional to observed codon counts
    @param subs_counts: observed substitution counts
    """
    score_matrix_transpose = P.T * v
    return algopy.sum(subs_counts * algopy.log(score_matrix_transpose))

def inner_eval_f(
        theta,
        pre_Q_postfix,
        subs_counts, v,
        ts, tv, syn, nonsyn,
        ):
    """
    This function is meant to be optimized with the help of algopy and ncg.
    @param theta: vector of unconstrained free variables
    @param pre_Q_postfix: this has estimates from the outer ML loop
    @param subs_counts: empirical substitution counts
    @param v: empirical codon distribution
    @param ts: precomputed nucleotide transition mask
    @param tv: precomputed nucleotide transversion mask
    @param syn: precomputed synonymous codon change mask
    @param nonsyn: precomputed non-synonymous codon change mask
    @return: negative log likelihood
    """
    log_mu = theta[0]
    log_kappa = theta[1]
    log_omega = theta[2]
    pre_Q_prefix = get_Q_prefix(
            ts, tv, syn, nonsyn,
            log_mu, log_kappa, log_omega)
    pre_Q = pre_Q_prefix * pre_Q_postfix
    Q = pre_Q - algopy.diag(algopy.sum(pre_Q, axis=1))
    P = algopy.expm(Q)
    return -get_log_likelihood(P, v, subs_counts)

def eval_f_unconstrained(
        theta,
        subs_counts, log_counts, v,
        h,
        ts, tv, syn, nonsyn, compo, asym_compo,
        ):
    """
    @param theta: vector of unconstrained free variables
    """
    d = theta[0]
    log_nt_weights = numpy.zeros(4, dtype=theta)
    log_nt_weights[0] = theta[1]
    log_nt_weights[1] = theta[2]
    log_nt_weights[2] = theta[3]
    log_nt_weights[3] = 0
    # laboriously construct this matrix
    pre_Q_postfix = get_Q_postfix_recessivity(
            compo, asym_compo,
            h,
            log_counts,
            d, log_nt_weights)
    #
    # construct the transition matrix
    Q = get_Q_unconstrained(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, d, log_nt_weights)
    P = algopy.expm(Q)
    #
    g = functools.partial(eval_grad, f)
    h = functools.partial(eval_hess, f)
    results = scipy.optimize.fmin_ncg(
            f,
            theta,
            args=fmin_args,
            fprime=g,
            fhess=h,
            maxiter=10000,
            avextol=1e-6,
            full_output=True,
            disp=True,
            retall=True,
            )
    xopt = results[0]
    yopt = results[1]
    print 'xopt:', xopt
    print 'yopt:', yopt
    return yopt


def submain_unconstrained_dominance(args):
    #
    # Precompute some ndarrays
    # according to properties of DNA and the genetic code.
    if args.mtdna or args.force_mtcode:
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
    for a, b in zip(codons, codon_counts):
        print a, ':', b
    print 'raw codon total:', numpy.sum(codon_counts)
    print 'raw codon counts:', codon_counts
    codon_counts = codon_counts[:len(codons)]
    print 'non-stop codon total:', numpy.sum(codon_counts)
    subs_counts = subs_counts[:len(codons), :len(codons)]
    v = codon_counts / float(numpy.sum(codon_counts))
    log_counts = numpy.log(codon_counts)
    #
    if args.disease == 'unconstrained':
        if args.integrate == 'quadrature':
            h = get_fixation_unconstrained_quad
        elif args.integrate == 'special':
            h = get_fixation_unconstrained
        else:
            raise Exception
    else:
        raise Exception
    #
    # predefine some plausible parameters but not the scaling parameter
    log_mu = 0
    log_kappa = 1
    log_omega = -3
    d = 0.5
    log_nt_weights = numpy.zeros(4)
    #
    # get the rate matrix associated with the initial guess
    Q = get_Q_unconstrained(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, d, log_nt_weights)
    #
    # get the minimum expected number of substitutions between codons
    mu_empirical = npcodon.get_lb_expected_subs(ham, subs_counts)
    mu_implied = -numpy.sum(numpy.diag(Q) * v)
    log_mu = math.log(mu_empirical) - math.log(mu_implied)
    print 'lower bound on expected mutations per codon site:', mu_empirical
    print
    # construct the initial guess
    theta = numpy.array([
        log_mu,
        log_kappa,
        log_omega,
        d,
        0,
        0,
        0,
        ])
    #
    # get the log likelihood associated with the initial guess
    fmin_args = (
            subs_counts, log_counts, v,
            h,
            ts, tv, syn, nonsyn, compo, asym_compo,
            )
    initial_cost = eval_f_unconstrained(theta, *fmin_args)
    print 'negative log likelihood of initial guess:',
    print initial_cost
    print
    print 'entropy bound on negative log likelihood:',
    print npcodon.get_lb_neg_ll(subs_counts)
    print
    f = eval_f_unconstrained
    results = scipy.optimize.fmin(
            f,
            theta,
            args=fmin_args,
            maxfun=10000,
            maxiter=10000,
            xtol=1e-8,
            ftol=1e-8,
            full_output=True,
            )

    print 'results:', results
    xopt = results[0]
    print 'optimal solution vector:', xopt
    print 'exp optimal solution vector:', numpy.exp(xopt)
    print
    print 'inverse of hessian:'
    print scipy.linalg.inv(h(xopt, *fmin_args))
    print


def eval_grad(f, theta, *args):
    """
    Compute the gradient of f in the forward mode of automatic differentiation.
    """
    theta = algopy.UTPM.init_jacobian(theta)
    retval = f(theta, *args)
    return algopy.UTPM.extract_jacobian(retval)

def eval_hess(f, theta, *args):
    """
    Compute the hessian of f in the forward mode of automatic differentiation.
    """
    theta = algopy.UTPM.init_hessian(theta)
    retval = f(theta, *args)
    return algopy.UTPM.extract_hessian(len(theta), retval)


def main(args):
    return submain_unconstrained_dominance(args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--mtdna',
            action='store_true',
            help='read the mtdna file from the website of Ziheng Yang')
    parser.add_argument(
            '--infile',
            help='codon alignment input file using the format of Ziheng Yang',
            required=True)
    parser.add_argument(
            '--t1',
            default='mm8',
            help='name of first taxon')
    parser.add_argument(
            '--t2',
            default='rn4',
            help='name of second taxon')
    main(parser.parse_args())

