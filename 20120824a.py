"""
Check continuous approximations of sojourn density.

Show the proportion of dimorphic time spent at each
allele frequency after a preferred
(or unpreferred if selection is negative)
mutation has occurred in a population fixed for the non-mutant allele
and before the mutant has gone to fixation or has been lost.
"""

from StringIO import StringIO
import time
import math

import numpy as np
from scipy import integrate
from scipy import special

import Form
import FormOut
import MatrixUtil
import StatsUtil
import kaizeng
import wrightfisher
import Util
import RUtil
from RUtil import mk_call_str
import wfengine

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size',
                5, low=3, high=40),
            Form.Float('gamma', 'scaled genic selection value', 1.5),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def sojourn_indefinite(x, g):
    """
    integral of 2*expm1(-a*(1-x)) / (expm1(-a)*x*(1-x))
    limit x -> 1 of (-2/expm1(a))*(-exp(a)*Ei(a*(x-1)) + Ei(a*x) + exp(a)*log((1-x)/x))
    """
    if g:
        prefix = -2 / math.expm1(g)
        suffix = 0
        suffix += -math.exp(g)*special.expi(g*(x-1))
        suffix += special.expi(g*x)
        suffix += math.exp(g)*math.log((1-x)/x)
        return prefix * suffix
    else:
        return 2*math.log(x)

def sojourn_definite(x0, x1, g):
    return sojourn_indefinite(x1, g) - sojourn_indefinite(x0, g)

def get_response_content(fs):
    N_diploid = fs.N_diploid
    N = N_diploid * 2
    k = 2
    gamma = fs.gamma
    # define the fitnesses and the selection value
    f0 = 1.0
    f1 = 1.0 - gamma / N
    s = 1 - f1 / f0
    if f1 <= 0:
        raise ValueError('the extreme selection caused a non-positive fitness')
    # get a wright fisher transition matrix
    P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
    """
    # condition on no fixation
    for i in range(N):
        P[i] /= 1 - P[i, N]
    # remove the fixed state from the transition matrix
    P = P[:N, :N]
    """
    # add mutations
    P[0, 0] = 0
    P[0, 1] = 1
    P[N, N] = 0
    P[N, 1] = 1
    # compute the stationary distribution
    v = MatrixUtil.get_stationary_distribution(P)
    # get the distribution over dimorphic states
    h = v[1:-1]
    h /= np.sum(h)
    # look at continuous approximations
    w = np.zeros(N+1)
    for i in range(1, N):
        """
        x = i / float(N)
        top = 2 * math.expm1(-gamma*(1-x))
        bottom = math.expm1(-gamma) * x * (1-x)
        value = top / bottom
        """
        #"""
        x0 = (i - 0.5) / float(N)
        #x0 = (i - 0.0) / float(N)
        x1 = (i + 0.5) / float(N)
        value = sojourn_definite(x0, x1, gamma)
        #"""
        w[i] = value
    w = w[1:-1]
    w /= np.sum(w)
    # get the array for the R plot
    arr = [h.tolist(), w.tolist()]
    # define the r script
    out = StringIO()
    print >> out, 'title.string <- "allele 1 vs allele 2"'
    print >> out, 'mdat <-', RUtil.matrix_to_R_string(arr)
    print >> out, mk_call_str(
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
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

