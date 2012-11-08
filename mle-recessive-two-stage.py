"""
This is intended to be a faster version of the similarly named script.

It is an exercise in premature optimization.
Or maybe it is not premature.
"""

import math
import argparse
import functools

import numpy
import scipy
import scipy.optimize
import scipy.linalg
import algopy

import jeffopt
import kimrecessive
import npcodon
import yangdata


##########################################################################
# These two functions could possibly be moved out of this module.
# But they could depend on algopy in other contexts.
# For example if you need the gradient with respect to the nucleotide weights,
# then these functions would need to be algopy-enabled.
# But in this script that is not the case.

def get_selection_F(log_counts, compo, log_nt_weights):
    """
    The F and S notation is from Yang and Nielsen 2008.
    Note that three of the four log nt weights are free parameters.
    One of the four log weights is zero and the other three
    are free parameters to be estimated jointly in the
    maximimum likelihood search,
    so this function is inside an optimization loop.
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
# In this script, these functions do not require algopy.
# Instead we use numerical quadrature and do not explicitly track
# the gradient with respect to the d parameter
# or with respect to any of the mutational process nucleotide
# equilibrium parameters that help to define S.

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

def get_Q_suffix(
        compo, asym_compo,
        log_counts,
        d, log_nt_weights):
    """
    This is specific to the model of recessivity.
    """
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F)
    pre_Q_suffix = numpy.exp(numpy.dot(asym_compo, log_nt_weights))
    pre_Q_suffix *= get_fixation_unconstrained_quad(S, d)
    return pre_Q_suffix


##########################################################################
# AlgoPy stuff involving parameters.
# The first few functions can be reused for all models
# regardless of how they treat mutational exchangeability.
# Some of the subsequent functions depend on how this is treated.


def get_Q(pre_Q_prefix, pre_Q_suffix):
    """
    @param pre_Q_prefix: this is an algopy aware ndarray
    @param pre_Q_suffix: this is a numpy ndarray
    @return: an algopy aware ndarray
    """
    pre_Q = pre_Q_prefix * pre_Q_suffix
    Q = pre_Q - algopy.diag(algopy.sum(pre_Q, axis=1))
    return Q

def help_log_lik(a, b):
    #return a*b
    return b*a

def get_log_likelihood(P, v, subs_counts):
    """
    The stationary distribution of P is empirically derived.
    It is proportional to the codon counts by construction.
    @param P: a transition matrix using codon counts and free parameters
    @param v: stationary distribution proportional to observed codon counts
    @param subs_counts: observed substitution counts
    """
    #FIXME: this function has been broken up for optimization debugging
    #score_matrix = P.T * v
    score_matrix = algopy.dot(algopy.diag(v), P)
    log_score_matrix = algopy.log(score_matrix)
    #log_likelihoods = subs_counts * log_score_matrix
    log_likelihoods = help_log_lik(subs_counts, log_score_matrix)
    log_likelihood = algopy.sum(log_likelihoods)
    return log_likelihood

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

def get_Q_prefix_gtr(
        gtr, syn, nonsyn,
        log_mu, log_gtr_exch, log_omega):
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
    gtr_exch = algopy.exp(log_gtr_exch)
    omega = algopy.exp(log_omega)
    return mu * algopy.dot(gtr, gtr_exch) * (omega * nonsyn + syn)

def inner_eval_f(
        theta,
        pre_Q_suffix,
        subs_counts, v,
        ts, tv, syn, nonsyn,
        ):
    """
    This function is meant to be optimized with the help of algopy and ncg.
    @param theta: vector of unconstrained free variables
    @param pre_Q_suffix: this has estimates from the outer ML loop
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
    Q = get_Q(pre_Q_prefix, pre_Q_suffix)
    P = algopy.expm(Q)
    return -get_log_likelihood(P, v, subs_counts)

def inner_eval_f_gtr(
        theta,
        pre_Q_suffix,
        subs_counts, v,
        gtr, syn, nonsyn,
        ):
    """
    This function is meant to be optimized with the help of algopy and ncg.
    @param theta: vector of unconstrained free variables
    @param pre_Q_suffix: this has estimates from the outer ML loop
    @param subs_counts: empirical substitution counts
    @param v: empirical codon distribution
    @param ts: precomputed nucleotide transition mask
    @param tv: precomputed nucleotide transversion mask
    @param syn: precomputed synonymous codon change mask
    @param nonsyn: precomputed non-synonymous codon change mask
    @return: negative log likelihood
    """
    log_mu = theta[0]
    log_gtr_exch = algopy.zeros(6, dtype=theta)
    log_gtr_exch[0] = theta[1]
    log_gtr_exch[1] = theta[2]
    log_gtr_exch[2] = theta[3]
    log_gtr_exch[3] = theta[4]
    log_gtr_exch[4] = theta[5]
    log_gtr_exch[5] = 0
    log_omega = theta[6]
    pre_Q_prefix = get_Q_prefix_gtr(
            gtr, syn, nonsyn,
            log_mu, log_gtr_exch, log_omega)
    Q = get_Q(pre_Q_prefix, pre_Q_suffix)
    P = algopy.expm(Q)
    return -get_log_likelihood(P, v, subs_counts)

def eval_f_unconstrained(
        theta,
        mu_empirical, subs_counts, log_counts, v,
        ts, tv, syn, nonsyn, compo, asym_compo,
        boxed_guess,
        ):
    """
    This function depends on the recessivity model.
    Nothing passed into this function is algopy aware.
    @param theta: vector of unconstrained free variables
    @return: negative log likelihood
    """
    print 'outer params:', theta
    print 'outer params exp:', numpy.exp(theta)
    d = theta[0]
    log_nt_weights = numpy.zeros(4)
    log_nt_weights[0] = theta[1]
    log_nt_weights[1] = theta[2]
    log_nt_weights[2] = theta[3]
    log_nt_weights[3] = 0
    # construct the suffix matrix using slow numerical integration
    pre_Q_suffix = get_Q_suffix(
            compo, asym_compo,
            log_counts,
            d, log_nt_weights,
            )
    #FIXME: use info from prev iterations to construct this guess
    # Construct an initial guess for the inner optimization.
    # log of generic scaling parameter
    # log of transition vs. transversion exchangeability ratio
    # log of nonsynonymous vs. synonymous exchangeability ratio
    if boxed_guess[0] is None:
        log_mu = 0.0
        log_kappa = 1.0
        log_omega = -1.0
        # re-estimate the generic scaling parameter
        pre_Q_prefix = get_Q_prefix(
                ts, tv, syn, nonsyn,
                log_mu, log_kappa, log_omega)
        Q = get_Q(pre_Q_prefix, pre_Q_suffix)
        mu_implied = -numpy.dot(numpy.diag(Q), v)
        log_mu = math.log(mu_empirical) - math.log(mu_implied)
        inner_guess = numpy.array([log_mu, log_kappa, log_omega])
    else:
        inner_guess = boxed_guess[0]
    # get conditional max likelihood estimates of the three inner parameters
    fmin_args = (
            pre_Q_suffix,
            subs_counts, v,
            ts, tv, syn, nonsyn,
            )
    f = inner_eval_f
    g = functools.partial(eval_grad, f)
    h = functools.partial(eval_hess, f)
    results = scipy.optimize.fmin_ncg(
            f,
            inner_guess,
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
    print 'inner xopt:', xopt
    print 'inner xopt exp:', numpy.exp(xopt)
    print 'inner neg log likelihood:', yopt
    boxed_guess[0] = xopt
    return yopt

def eval_f_unconstrained_gtr(
        theta,
        mu_empirical, subs_counts, log_counts, v,
        gtr, syn, nonsyn, compo, asym_compo,
        boxed_guess,
        ):
    """
    This function depends on the recessivity model.
    Nothing passed into this function is algopy aware.
    @param theta: vector of unconstrained free variables
    @return: negative log likelihood
    """
    print 'outer params:', theta
    print 'outer params exp:', numpy.exp(theta)
    d = theta[0]
    log_nt_weights = numpy.zeros(4)
    log_nt_weights[0] = theta[1]
    log_nt_weights[1] = theta[2]
    log_nt_weights[2] = theta[3]
    log_nt_weights[3] = 0
    # construct the suffix matrix using slow numerical integration
    pre_Q_suffix = get_Q_suffix(
            compo, asym_compo,
            log_counts,
            d, log_nt_weights,
            )
    # Construct an initial guess for the inner optimization.
    if boxed_guess[0] is None:
        log_mu = 0.0
        log_gtr_exch = numpy.zeros(6)
        log_omega = -1.0
        # re-estimate the generic scaling parameter
        pre_Q_prefix = get_Q_prefix_gtr(
                gtr, syn, nonsyn,
                log_mu, log_gtr_exch, log_omega)
        Q = get_Q(pre_Q_prefix, pre_Q_suffix)
        mu_implied = -numpy.dot(numpy.diag(Q), v)
        log_mu = math.log(mu_empirical) - math.log(mu_implied)
        inner_guess = numpy.array([
            log_mu,
            0, 0, 0, 0, 0,
            #log_gtr_exch[0],
            #log_gtr_exch[1],
            #log_gtr_exch[2],
            #log_gtr_exch[3],
            #log_gtr_exch[4],
            log_omega,
            ])
    else:
        inner_guess = boxed_guess[0]
    # get conditional max likelihood estimates of the three inner parameters
    fmin_args = (
            pre_Q_suffix,
            subs_counts, v,
            gtr, syn, nonsyn,
            )
    f = inner_eval_f_gtr
    g = functools.partial(eval_grad, f)
    h = functools.partial(eval_hess, f)
    results = scipy.optimize.fmin_ncg(
            f,
            inner_guess,
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
    print 'inner xopt:', xopt
    print 'inner xopt exp:', numpy.exp(xopt)
    print 'inner neg log likelihood:', yopt
    boxed_guess[0] = xopt
    return yopt


def submain_unconstrained_dominance(args):
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
    # get the minimum expected number of substitutions between codons
    mu_empirical = npcodon.get_lb_expected_subs(ham, subs_counts)
    print 'lower bound on expected mutations per codon site:', mu_empirical
    print
    print 'entropy lower bound on negative log likelihood:',
    print npcodon.get_lb_neg_ll(subs_counts)
    print
    #
    # initialize parameter value guesses
    d = 0.5
    theta = numpy.array([
        d,
        0, 0, 0,
        ], dtype=float)
    boxed_guess = [None]
    fmin_args = (
            mu_empirical, subs_counts, log_counts, v,
            ts, tv, syn, nonsyn, compo, asym_compo,
            boxed_guess,
            )
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
    #print 'inverse of hessian:'
    #print scipy.linalg.inv(h(xopt, *fmin_args))
    #print

def submain_unconstrained_dominance_gtr(args):
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
    # get the minimum expected number of substitutions between codons
    mu_empirical = npcodon.get_lb_expected_subs(ham, subs_counts)
    print 'lower bound on expected mutations per codon site:', mu_empirical
    print
    print 'entropy lower bound on negative log likelihood:',
    print npcodon.get_lb_neg_ll(subs_counts)
    print
    #
    # initialize parameter value guesses
    d = 0.5
    theta = numpy.array([
        d,
        0, 0, 0,
        ], dtype=float)
    boxed_guess = [None]
    fmin_args = (
            mu_empirical, subs_counts, log_counts, v,
            gtr, syn, nonsyn, compo, asym_compo,
            boxed_guess,
            )
    f = eval_f_unconstrained_gtr
    """
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
    """
    results = scipy.optimize.minimize(
            f,
            theta,
            args=fmin_args,
            method='Nelder-Mead',
            )
    print 'results:', results
    xopt = results[0]
    print 'optimal solution vector:', xopt
    print 'exp optimal solution vector:', numpy.exp(xopt)
    print
    #print 'inverse of hessian:'
    #print scipy.linalg.inv(h(xopt, *fmin_args))
    #print


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
    if args.mutexch == 'hky':
        return submain_unconstrained_dominance(args)
    elif args.mutexch == 'gtr':
        return submain_unconstrained_dominance_gtr(args)
    else:
        raise Exception


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--mutexch',
            choices=('hky', 'gtr'),
            default='hky',
            help='model of mutational exchangeability')
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

