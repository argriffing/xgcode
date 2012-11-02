"""
This is currently a script but later it should be web accessible.
First version is October 4 2012.
It is for testing a pure recessive disease codon preference model
against the genic model of Yang and Nielsen 2008.
This script does not use so many biological libraries
because for some reason I am adopting a more matlab-like style.
All matrices in this script are ndarrays as opposed to actual python matrices.
.
Added October 30 2012 --
This script has been adapted as the codon_model.py example in the python
automatic differentiation package called algopy.
Now I am backporting those algopy adaptations back into this more complicated
mle-recessive.py script to make use of the algopy features.
This means that some of the arrays that were formerly numpy ndarrays
are now algopy UTPM objects that mimic arrays and that carry
nonzero-degree Taylor terms through the log likelihood computation.
These higher order Taylor terms have a couple of uses.
One of the uses is to provide gradient and hessian information to
relatively sophisticated likelihood maximization algorithms
like the scipy truncated newton or conjugate gradient methods.
Another use of these higher order terms is to provide estimates of
parameter uncertainty by taking the inverse of the hessian of the
negative log likelihood at the max likelihood point as the covariance matrix
that approximates something like a multivariate normal
joint posterior distribution of the parameters.
Although this interpretation as a posterior distribution is probably
horribly flawed unless it assumes an implicit prior distribution that
trivially forces this to be a posterior distribution.
"""

import string
import math
import argparse
from itertools import product

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
# these two functions are pure numpy and speed does not matter


def get_lb_neg_ll(subs_counts):
    """
    Get the lower bound negative log likelihood.
    It uses an entropy-based calculation.
    @param subs_counts: codon substitution counts
    """
    nstates = subs_counts.shape[0]
    counts = []
    for i in range(nstates):
        for j in range(nstates):
            if i < j:
                c = subs_counts[i, j] + subs_counts[j, i]
                counts.append(c)
            elif i == j:
                c = subs_counts[i, j]
                counts.append(c)
    # return the entropy of the unordered pair count vector
    probs = numpy.array(counts, dtype=float) / numpy.sum(counts)
    return numpy.sum(-c*math.log(p) for c, p in zip(counts, probs) if c)

def get_lb_expected_subs(ham, subs_counts):
    """
    Get the lower bound of expected substitutions.
    This is expected minimum possible number of substitutions
    between aligned codons.
    """
    return numpy.sum(ham * subs_counts) / float(numpy.sum(subs_counts))


##########################################################################
# algopy stuff involving parameters

def get_fixation_knudsen(S):
    """
    This is +gwF = 1/2.
    """
    return 1. / kimrecessive.denom_knudsen(0.5*S)

def get_fixation_genic(S):
    return 1. / kimrecessive.denom_genic_a(0.5*S)

def get_fixation_recessive_disease(S):
    sign_S = algopy.sign(S)
    H = algopy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_piecewise(
                    0.5*S[i, j], sign_S[i, j])
    return H

def get_fixation_dominant_disease(S):
    sign_S = algopy.sign(S)
    H = algopy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_piecewise(
                    0.5*S[i, j], -sign_S[i, j])
    return H

def get_fixation_unconstrained(S, d):
    sign_S = algopy.sign(S)
    D = d * sign_S
    H = algopy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_piecewise(
                    0.5*S[i, j], D[i, j])
    return H

def get_fixation_unconstrained_kb(S, d, log_kb):
    """
    This uses the Kacser and Burns effect instead of the sign function.
    """
    soft_sign_S = algopy.tanh(algopy.exp(log_kb)*S)
    D = d * soft_sign_S
    H = algopy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_piecewise(
                    0.5*S[i, j], D[i, j])
    return H

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
    return log_counts - algopy.dot(compo, log_nt_weights)

def get_selection_S(F):
    """
    The F and S notation is from Yang and Nielsen 2008.
    @param F: a selection value for each codon, up to an additive constant
    @return: selection differences F_j - F_i, also known as S_ij
    """
    e = algopy.ones_like(F)
    return algopy.outer(e, F) - algopy.outer(F, e)

def get_Q(
        ts, tv, syn, nonsyn, compo, asym_compo,
        h,
        log_counts,
        log_mu, log_kappa, log_omega, log_nt_weights):
    """
    Notation is from Yang and Nielsen 2008.
    The first group of args consists of precomputed ndarrays.
    The second group is only the fixation function.
    The third group consists of empirically (non-free) estimated parameters.
    The fourth group depends only on free parameters.
    @param ts: indicator for transition
    @param tv: indicator for transversion
    @param syn: indicator for synonymous codons
    @param nonsyn: indicator for nonsynonymous codons
    @param compo: site independent nucleotide composition per codon
    @param asym_compo: tensor from get_asym_compo function
    @param h: fixation function
    @param log_counts: empirically counted codons in the data set
    @param log_mu: free param for scaling
    @param log_kappa: free param for transition transversion rate distinction
    @param log_omega: free param for syn nonsyn rate distinction
    @param log_nt_weights: mostly free param array for mutation equilibrium
    @return: rate matrix
    """
    mu = algopy.exp(log_mu)
    kappa = algopy.exp(log_kappa)
    omega = algopy.exp(log_omega)
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F)
    pre_Q = mu * (kappa * ts + tv) * (omega * nonsyn + syn) * algopy.exp(
            algopy.dot(asym_compo, log_nt_weights)) * h(S)
    Q = pre_Q - algopy.diag(algopy.sum(pre_Q, axis=1))
    return Q

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

def get_Q_unconstrained_kb(
        ts, tv, syn, nonsyn, compo, asym_compo,
        h,
        log_counts,
        log_mu, log_kappa, log_omega, d, log_kb, log_nt_weights):
    """
    This adds yet another parameter.
    """
    mu = algopy.exp(log_mu)
    kappa = algopy.exp(log_kappa)
    omega = algopy.exp(log_omega)
    F = get_selection_F(log_counts, compo, log_nt_weights)
    S = get_selection_S(F)
    pre_Q = mu * (kappa * ts + tv) * (omega * nonsyn + syn) * algopy.exp(
            algopy.dot(asym_compo, log_nt_weights)) * h(S, d, log_kb)
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
    return algopy.sum(subs_counts * algopy.log(P.T * v))

def eval_f(
        theta,
        subs_counts, log_counts, v,
        h,
        ts, tv, syn, nonsyn, compo, asym_compo,
        ):
    """
    The function formerly known as minimize-me.
    @param theta: length six unconstrained vector of free variables
    """
    # unpack theta
    log_mu = theta[0]
    log_kappa = theta[1]
    log_omega = theta[2]
    log_nt_weights = algopy.zeros(4, dtype=theta)
    log_nt_weights[0] = theta[3]
    log_nt_weights[1] = theta[4]
    log_nt_weights[2] = theta[5]
    log_nt_weights[3] = 0
    #
    # construct the transition matrix
    Q = get_Q(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, log_nt_weights)
    P = algopy.expm(Q)
    #
    # return the neg log likelihood
    return -get_log_likelihood(P, v, subs_counts)

def eval_f_unconstrained(
        theta,
        subs_counts, log_counts, v,
        h,
        ts, tv, syn, nonsyn, compo, asym_compo,
        ):
    """
    No dominance/recessivity constraint.
    @param theta: length seven unconstrained vector of free variables
    """
    # unpack theta
    log_mu = theta[0]
    log_kappa = theta[1]
    log_omega = theta[2]
    d = theta[3]
    log_nt_weights = algopy.zeros(4, dtype=theta)
    log_nt_weights[0] = theta[4]
    log_nt_weights[1] = theta[5]
    log_nt_weights[2] = theta[6]
    log_nt_weights[3] = 0
    #
    # construct the transition matrix
    Q = get_Q_unconstrained(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, d, log_nt_weights)
    P = algopy.expm(Q)
    #
    # return the neg log likelihood
    return -get_log_likelihood(P, v, subs_counts)

def eval_f_unconstrained_kb(
        theta,
        subs_counts, log_counts, v,
        h,
        ts, tv, syn, nonsyn, compo, asym_compo,
        ):
    """
    No dominance/recessivity constraint.
    @param theta: length seven unconstrained vector of free variables
    """
    # unpack theta
    log_mu = theta[0]
    log_kappa = theta[1]
    log_omega = theta[2]
    d = theta[3]
    log_kb = theta[4]
    log_nt_weights = algopy.zeros(4, dtype=theta)
    log_nt_weights[0] = theta[5]
    log_nt_weights[1] = theta[6]
    log_nt_weights[2] = theta[7]
    log_nt_weights[3] = 0
    #
    # construct the transition matrix
    Q = get_Q_unconstrained_kb(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, d, log_kb, log_nt_weights)
    P = algopy.expm(Q)
    #
    # return the neg log likelihood
    return -get_log_likelihood(P, v, subs_counts)

def eval_grad_f(theta, *args):
    """
    compute the gradient of f in the forward mode of AD
    """
    theta = algopy.UTPM.init_jacobian(theta)
    retval = eval_f(theta, *args)
    return algopy.UTPM.extract_jacobian(retval)

def eval_hess_f(theta, *args):
    """
    compute the hessian of f in the forward mode of AD
    """
    theta = algopy.UTPM.init_hessian(theta)
    retval = eval_f(theta, *args)
    return algopy.UTPM.extract_hessian(len(theta), retval)

def eval_grad_f_unconstrained(theta, *args):
    """
    compute the gradient of f in the forward mode of AD
    No dominance/recessivity constraint.
    """
    theta = algopy.UTPM.init_jacobian(theta)
    retval = eval_f_unconstrained(theta, *args)
    return algopy.UTPM.extract_jacobian(retval)

def eval_hess_f_unconstrained(theta, *args):
    """
    compute the hessian of f in the forward mode of AD
    No dominance/recessivity constraint.
    """
    theta = algopy.UTPM.init_hessian(theta)
    retval = eval_f_unconstrained(theta, *args)
    return algopy.UTPM.extract_hessian(len(theta), retval)

def eval_grad_f_unconstrained_kb(theta, *args):
    """
    compute the gradient of f in the forward mode of AD
    No dominance/recessivity constraint.
    """
    theta = algopy.UTPM.init_jacobian(theta)
    retval = eval_f_unconstrained_kb(theta, *args)
    return algopy.UTPM.extract_jacobian(retval)

def eval_hess_f_unconstrained_kb(theta, *args):
    """
    compute the hessian of f in the forward mode of AD
    No dominance/recessivity constraint.
    """
    theta = algopy.UTPM.init_hessian(theta)
    retval = eval_f_unconstrained_kb(theta, *args)
    return algopy.UTPM.extract_hessian(len(theta), retval)

def submain_unconstrained_dominance_kb(args):
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
    codon_counts = np.sum(subs_counts, axis=0) + np.sum(subs_counts, axis=1)
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
    if args.disease == 'kacser':
        h = get_fixation_unconstrained_kb
    else:
        raise Exception
    #
    # predefine some plausible parameters but not the scaling parameter
    log_mu = 0
    log_kappa = 1
    log_omega = -3
    d = 0.5
    log_kb = 0
    log_nt_weights = numpy.zeros(4)
    #
    # get the rate matrix associated with the initial guess
    Q = get_Q_unconstrained_kb(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, d, log_kb, log_nt_weights)
    #
    # get the minimum expected number of substitutions between codons
    mu_empirical = get_lb_expected_subs(ham, subs_counts)
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
        log_kb,
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
    initial_cost = eval_f_unconstrained_kb(theta, *fmin_args)
    print 'negative log likelihood of initial guess:',
    print initial_cost
    print
    print 'entropy bound on negative log likelihood:',
    print get_lb_neg_ll(subs_counts)
    print
    #
    # search for the minimum negative log likelihood over multiple parameters
    if args.fmin == 'ncg':
        result = scipy.optimize.fmin_ncg(
                eval_f_unconstrained_kb,
                theta,
                args=fmin_args,
                fprime=eval_grad_f_unconstrained_kb,
                fhess=eval_hess_f_unconstrained_kb,
                maxiter=10000,
                full_output=True,
                disp=True,
                retall=True,
                )
    else:
        raise Exception
    print 'results:', result
    xopt = result[0]
    print 'optimal solution vector:', xopt
    print 'exp optimal solution vector:', numpy.exp(xopt)
    print
    print 'inverse of hessian:'
    print scipy.linalg.inv(eval_hess_f_unconstrained_kb(xopt, *fmin_args))
    print


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
    codon_counts = np.sum(subs_counts, axis=0) + np.sum(subs_counts, axis=1)
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
        h = get_fixation_unconstrained
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
    mu_empirical = get_lb_expected_subs(ham, subs_counts)
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
    print get_lb_neg_ll(subs_counts)
    print
    #
    # search for the minimum negative log likelihood over multiple parameters
    if args.fmin == 'ncg':
        result = scipy.optimize.fmin_ncg(
                eval_f_unconstrained,
                theta,
                args=fmin_args,
                fprime=eval_grad_f_unconstrained,
                fhess=eval_hess_f_unconstrained,
                maxiter=10000,
                full_output=True,
                disp=True,
                retall=True,
                )
    else:
        raise Exception
    print 'results:', result
    xopt = result[0]
    print 'optimal solution vector:', xopt
    print 'exp optimal solution vector:', numpy.exp(xopt)
    print
    print 'inverse of hessian:'
    print scipy.linalg.inv(eval_hess_f_unconstrained(xopt, *fmin_args))
    print



def submain_constrained_dominance(args):
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
    codon_counts = numpy.sum(subs_counts, axis=0) + numpy.sum(
            subs_counts, axis=1)
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
    if args.disease == 'genic':
        h = get_fixation_genic
    elif args.disease == 'recessive':
        h = get_fixation_recessive_disease
    elif args.disease == 'dominant':
        h = get_fixation_dominant_disease
    else:
        raise Exception
    #
    # predefine some plausible parameters but not the scaling parameter
    log_mu = 0
    log_kappa = 1
    log_omega = -3
    log_nt_weights = numpy.zeros(4)
    #
    # get the rate matrix associated with the initial guess
    Q = get_Q(
            ts, tv, syn, nonsyn, compo, asym_compo,
            h,
            log_counts,
            log_mu, log_kappa, log_omega, log_nt_weights)
    #
    # get the minimum expected number of substitutions between codons
    mu_empirical = get_lb_expected_subs(ham, subs_counts)
    mu_implied = -numpy.sum(numpy.diag(Q) * v)
    log_mu = math.log(mu_empirical) - math.log(mu_implied)
    print 'lower bound on expected mutations per codon site:', mu_empirical
    print
    # construct the initial guess
    theta = numpy.array([
        log_mu,
        log_kappa,
        log_omega,
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
    initial_cost = eval_f(theta, *fmin_args)
    print 'negative log likelihood of initial guess:',
    print initial_cost
    print
    print 'entropy bound on negative log likelihood:',
    print get_lb_neg_ll(subs_counts)
    print
    #
    # search for the minimum negative log likelihood over multiple parameters
    if args.fmin == 'simplex':
        results = scipy.optimize.fmin(
                eval_f,
                theta,
                args=fmin_args,
                maxfun=10000,
                maxiter=10000,
                xtol=1e-8,
                ftol=1e-8,
                full_output=True,
                )
    elif args.fmin == 'bfgs':
        results = scipy.optimize.fmin_bfgs(
                eval_f,
                theta,
                args=fmin_args,
                maxiter=10000,
                full_output=True,
                )
    elif args.fmin == 'jeffopt':
        results = jeffopt.fmin_jeff_unconstrained(
                eval_f,
                theta,
                args=fmin_args,
                )
    elif args.fmin == 'ncg':
        results = scipy.optimize.fmin_ncg(
                eval_f,
                theta,
                args=fmin_args,
                fprime=eval_grad_f,
                fhess=eval_hess_f,
                maxiter=10000,
                full_output=True,
                disp=True,
                retall=True,
                )
    else:
        raise Exception
    print 'results:', results
    xopt = results[0]
    print 'optimal solution vector:', xopt
    print 'exp optimal solution vector:', numpy.exp(xopt)
    print


def main(args):
    #
    # Check reversibility of h functions with respect to F,
    # in the notation of Yang and Nielsen 2008.
    for h in (
            get_fixation_genic,
            get_fixation_recessive_disease,
            get_fixation_dominant_disease,
            get_fixation_knudsen,
            ):
        F = numpy.array([1.2, 2.3, 0, -1.1])
        S = get_selection_S(F)
        fixation = h(S)
        log_ratio = numpy.log(fixation / fixation.T)
        if not numpy.allclose(S, log_ratio):
            raise Exception((S, log_ratio))
    #
    # Do the main analysis.
    if args.disease == 'unconstrained':
        return submain_unconstrained_dominance(args)
    elif args.disease == 'kacser':
        return submain_unconstrained_dominance_kb(args)
    else:
        return submain_constrained_dominance(args)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
            '--fmin',
            choices=('simplex', 'bfgs', 'jeffopt', 'ncg'),
            default='simplex',
            help='nonlinear multivariate optimization')
    parser.add_argument(
            '--disease',
            choices=(
                'genic', 'recessive', 'dominant',
                'unconstrained', 'kacser'),
            default='genic',
            help='the mode of natural selection on unpreferred codons')
    parser.add_argument(
            '--mtdna',
            action='store_true',
            help='read the mtdna file from the website of Ziheng Yang')
    parser.add_argument(
            '--force-mtcode',
            action='store_true',
            help='use the mitochondrial genetic code and stop codons')
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

