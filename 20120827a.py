"""
Check Kimura approximation of expected fixation time given fixation.

This assumes a neutral allele.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import linalg

import Form
import FormOut
import MatrixUtil
import StatsUtil
import Util
import RUtil
from RUtil import mk_call_str
import wfengine
import kimura

def get_form():
    """
    @return: the body of a form
    """
    return [
            #Form.Integer('N_diploid', 'diploid population size',
                #5, low=3, high=40),
            Form.ImageFormat(),
            ]

def get_t1_exact(N_diploid):
    """
    @param N_diploid: number of diploid individuals in the population
    @return: a sequence of expected fixation times
    """
    s = 0
    N_haploid = N_diploid * 2
    P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
    # Create a transition matrix modified so that it is conditional
    # on eventual fixation.
    x = kimura.get_fixation_probabilities(P)
    A = (P * x)[1:, 1:]
    for i in range(N_haploid):
        A[i] /= np.sum(A[i])
    A[-1, -1] = 0
    A[-1, 0] = 1
    # Now use this conditional transition matrix
    # to set up a system of equations that will give
    # the expected number of generations until fixation
    # conditioned upon eventual fixation
    B = A - np.eye(N_haploid)
    b = -np.ones(N_haploid)
    B[-1] = np.zeros(N_haploid)
    B[-1, -1] = 1
    b[-1] = 0
    y = linalg.solve(B, b)
    return y.tolist()

def get_t1_approx(N_diploid):
    N_haploid = N_diploid * 2
    arr = []
    for i in range(1, N_haploid + 1):
        p = i / float(N_haploid)
        # Get the expected time until fixation given eventual fixation
        # when the population begins with proportion p mutants.
        if p == 1:
            t1 = 0
        else:
            t1 = -(1/p)*2*N_haploid*(1-p)*math.log(1-p)
        arr.append(t1)
    return arr

def get_form_out():
    return FormOut.Image('plot')

def get_response_content(fs):
    #
    N_diploid = 10
    N_haploid = 2 * N_diploid
    t1_20_exact = get_t1_exact(N_diploid)
    t1_20_approx = get_t1_approx(N_diploid)
    proportions_20 = np.linspace(0.0, 1.0, 2 * N_diploid + 1)[1:]
    #
    N_diploid = 5
    N_haploid = 2 * N_diploid
    t1_10_exact = get_t1_exact(N_diploid)
    t1_10_approx = get_t1_approx(N_diploid)
    proportions_10 = np.linspace(0.0, 1.0, 2 * N_diploid + 1)[1:]
    #
    # define the r script
    out = StringIO()
    print >> out, 'title.string <- "my title"'
    print >> out, 't1.20.approx <- c', str(tuple(t1_20_approx))
    print >> out, 't1.20.exact <- c', str(tuple(t1_20_exact))
    print >> out, 'p.20 <- c', str(tuple(proportions_20.tolist()))
    print >> out, mk_call_str('plot', 'p.20', 't1.20.approx',
            col='"red"', type='"l"', ylim='c(0,45)', xaxp='c(0,1,10)')
    print >> out, mk_call_str('points', 'p.20', 't1.20.exact',
            col='"blue"')
    #
    print >> out, 't1.10.approx <- c', str(tuple(t1_10_approx))
    print >> out, 't1.10.exact <- c', str(tuple(t1_10_exact))
    print >> out, 'p.10 <- c', str(tuple(proportions_10.tolist()))
    print >> out, mk_call_str('lines', 'p.10', 't1.10.approx',
            col='"red"')
    print >> out, mk_call_str('points', 'p.10', 't1.10.exact',
            col='"blue"')
    #
    """
            'barplot',
            'mdat',
            'legend.text=' + mk_call_str(
                'c',
                '"exact discrete distribution"',
                '"continuous approximation"',
                #'"two-allele large N limit"',
                #'"two-allele"',
                #'"four-allele without mutational bias"',
                #'"four-allele with mutational bias (kappa_{1,2}=2)"',
                #'"four-allele with mutational bias, large N limit"',
                ),
            'args.legend = list(x="topright", bty="n")',
            'names.arg = 1:%s' % (N-1),
            main='title.string',
            xlab='"frequency of allele 1"',
            ylab='"frequency"',
            col=mk_call_str(
                'c',
                #'"red"',
                #'"white"',
                '"black"',
                #'"gray"',
                '"red"',
                ),
            beside='TRUE',
            border='NA',
            )
    #print >> out, 'box()'
    """
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

