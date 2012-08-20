"""
Better approximate a figure from a publication by Kai Zeng 2010.

This approximation samples 10 alleles from a
fixed size population of diploid individuals.
The figure (fig. 2) in the publication may be for a large population limit.
"""

from StringIO import StringIO
import time
import math

import numpy as np

import Form
import FormOut
import MatrixUtil
import StatsUtil
import kaizeng
import wrightfisher
import Util
import RUtil
from RUtil import mk_call_str

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_big_diploid', 'diploid population size',
                100, low=5, high=150),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def get_two_allele_distribution(N_big, N_small, f0, f1):
    """
    Assumes small genic selection.
    Assumes small mutation.
    The mutational bias does not affect the distribution.
    @param N_big: total number of alleles in the population
    @param N_small: number of alleles sampled from the population
    @param f0: fitness of allele 0
    @param f1: fitness of allele 1
    @return: distribution over all non-fixed population states
    """
    # construct a transition matrix
    nstates = N_big + 1
    P = np.zeros((nstates, nstates))
    for i in range(nstates):
        p0, p1 = wrightfisher.genic_diallelic(f0, f1, i, N_big - i)
        if i == 0:
            P[i, 1] = 1.0
        elif i == N_big:
            P[i, N_big - 1] = 1.0
        else:
            for j in range(nstates):
                logp = StatsUtil.binomial_log_pmf(j, N_big, p0)
                P[i, j] = math.exp(logp)
    # find the stationary distribution
    v = MatrixUtil.get_stationary_distribution(P)
    # debugging
    if not np.allclose(np.sum(v), 1):
        raise ValueError
    if not all(x >= 0 for x in v):
        raise ValueError
    if not np.allclose(v, np.dot(v, P)):
        raise ValueError
    # return the stationary distribution conditional on dimorphism
    accum = np.zeros(N_small + 1)
    denominator = Util.choose(N_big, N_small)
    for i, p_state in enumerate(v):
        for j in range(N_small + 1):
            numerator = Util.choose(i, j) * Util.choose(N_big - i, N_small - j)
            if not numerator:
                continue
            accum[j] += (p_state * numerator) / denominator
    return accum[1:-1] / np.sum(accum[1:-1])

def get_response_content(fs):
    N_small = 10
    N_big = 2*fs.N_big_diploid
    if N_big < N_small:
        raise ValueError('use a larger diploid population size')
    denominator = Util.choose(N_big, N_small)
    k = 4
    gamma = 1.5
    params_list = [
            (0.008, 1, 1, 0, gamma, 1),
            (0.008, 2, 1, 0, gamma, 1)]
    allele_histograms = np.zeros((2, N_small + 1))
    for i, params in enumerate(params_list):
        mutation, selection = kaizeng.params_to_mutation_fitness(N_big, params)
        P = kaizeng.get_transition_matrix(N_big, k, mutation, selection)
        v = MatrixUtil.get_stationary_distribution(P)
        for state_index, counts in enumerate(kaizeng.gen_states(N_big, k)):
            if counts[0] and counts[1]:
                p_state = v[state_index]
                for j in range(1, N_small):
                    numerator = Util.choose(counts[0], j)
                    if not numerator: continue
                    numerator *= Util.choose(counts[1], N_small - j)
                    if not numerator: continue
                    x = (p_state * numerator) / denominator
                    allele_histograms[i, j] += x
    # Define the r table.
    # There are nine columns each corresponding to an allele frequency.
    # There are three rows each corresponding to a configuration.
    arr = []
    # Use the two allele approximation
    # from mcvean and charlesworth 1999 referred to by zeng 2011.
    # I'm not sure if I am using the right equation.
    gamma_0 = 0
    gamma_1 = 1.5
    s_0 = -gamma_0 / float(N_big)
    s_1 = -gamma_1 / float(N_big)
    hist = np.zeros(N_small+1)
    for i in range(1, N_small):
        x = i / float(N_small)
        hist[i] = math.exp(1*N_big*(s_0 - s_1)*x) / (x*(1-x))
    h = hist[1:-1]
    h /= np.sum(h)
    arr.append(h.tolist())
    # Use the exact two allele distribution.
    # Well, it is exact if I understand the right scaling
    # of the population size and fitnesses.
    f0 = 1.0
    f1 = 1.0 - gamma / N_big
    #f0 = 1.0 + gamma / N
    #f1 = 1.0
    #f0 = 1.0 + 1.5 / (4*N)
    #f1 = 1.0 - 1.5 / (4*N)
    h = get_two_allele_distribution(N_big, N_small, f0, f1)
    arr.append(h.tolist())
    # Get frequencies for the other two configurations
    for hist in allele_histograms:
        h = hist[1:-1]
        h /= np.sum(h)
        arr.append(h.tolist())
    # define the r script
    out = StringIO()
    print >> out, 'title.string <- "allele 1 vs allele 2, gamma = 1.5"'
    print >> out, 'mdat <-', RUtil.matrix_to_R_string(arr)
    print >> out, mk_call_str(
            'barplot',
            'mdat',
            'legend.text=' + mk_call_str(
                'c',
                '"two-allele diffusion limit"',
                '"two-allele"',
                '"four-allele without mutational bias"',
                '"four-allele with mutational bias kappa_{1,2}=2"',
                ),
            'args.legend = list(x="topleft", bty="n")',
            'names.arg = c(1,2,3,4,5,6,7,8,9)',
            main='title.string',
            xlab='"frequency of allele 1"',
            ylab='"frequency"',
            col=mk_call_str(
                'c',
                '"red"',
                '"white"',
                '"black"',
                '"gray"',
                ),
            beside='TRUE',
            )
    #print >> out, 'box()'
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data
