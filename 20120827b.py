"""
Check Kimura approximation of expected fixation time given fixation p0=1/10.

The initial frequency is assumed to be 1/10.
Plot over a range of selection strengths.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import linalg
from scipy import integrate

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

def get_t1_exact(N_mutants, N_diploid, s):
    """
    """
    N_haploid = N_diploid * 2
    p = N_mutants / float(N_haploid)
    #P = np.exp(wfengine.create_genic_diallelic(N_diploid, s))
    P = np.exp(wfengine.create_genic_diallelic_ohta(N_diploid, s))
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
    return y[1]

def f_psi(x, N_diploid, s):
    """
    The limits at s=0 and x=0 and x=1 are not implemented yet.
    Presumably this function will later be simplified,
    and then ultimately be used as a term in the integrand
    of an incomplete integral whose solution
    will possibly be a special function.
    """
    # integral of G from 0 to 1
    integral_of_G = -math.expm1(-4*N_diploid*s) / (4*N_diploid*s)
    # variance of allele frequency change per generation
    Vx = (x * (1-x)) / (2*N_diploid)
    # the G function evaluated at x
    Gx = math.exp(-4*N_diploid*s*x)
    # return the numerical value
    return (2 * integral_of_G) / (Vx * Gx)

def f_fix(x, N_diploid, s):
    """
    The limits at s=0 and x=0 and x=1 are not implemented yet.
    Ultimately this function will be used as a term of an integrand
    whose integral will possibly be a well known special function.
    """
    return math.expm1(-4*N_diploid*s*x) / math.expm1(-4*N_diploid*s)

def integrand_a(x, N_diploid, s):
    psi = f_psi(x, N_diploid, s)
    u = f_fix(x, N_diploid, s)
    return psi * u * (1-u)

def integrand_b(x, N_diploid, s):
    psi = f_psi(x, N_diploid, s)
    u = f_fix(x, N_diploid, s)
    return psi * u * u

def J_fix(x, S):
    """
    This is part of equation (17) of Kimura and Ohta.
    @param S: N_e * s
    """
    return math.expm1(-2*S*x) / math.expm1(-2*S)

def J1_integrand(x, S):
    """
    This is part of equation (17) of Kimura and Ohta.
    """
    a = math.expm1(2*S*x)
    b = math.exp(-2*S*x) - math.exp(-2*S)
    c = x * (1 - x)
    return (a * b) / c

def J2_integrand(x, S):
    """
    This is part of equation (17) of Kimura and Ohta.
    """
    a = math.expm1(2*S*x)
    b = math.expm1(-2*S*x)
    c = x / (1 - x)
    return -(a * b) / c

def get_t1_approx(p, N_diploid, s):
    """
    For now this will use numerical integration.
    """
    S = N_diploid * s
    #a = integrate.quad(J1_integrand, p, 1, args=(S,))[0]
    a0 = kimura.J1_indefinite_integral(p, 2*S)
    a1 = kimura.J1_indefinite_integral(1, 2*S)
    b = integrate.quad(J2_integrand, 0, p, args=(S,))[0]
    coeff = -2 / (s * math.expm1(-2*S))
    u = J_fix(p, S)
    return coeff * ((a1 - a0) + (u/(1-u)) * b)

def get_form_out():
    return FormOut.Image('plot')

def get_response_content(fs):
    N_diploid = 10
    N_mutants = 1
    N_haploid = N_diploid * 2
    p = N_mutants / float(N_haploid)
    t1_exact = []
    t1_approx = []
    Nes_values = range(1, 9)
    for Nes in Nes_values:
        s = Nes / float(N_diploid)
        t1_exact.append(get_t1_exact(N_mutants, N_diploid, s))
        t1_approx.append(get_t1_approx(p, N_diploid, s))
    # define the r script
    out = StringIO()
    print >> out, 'title.string <- "my title"'
    print >> out, 'Nes <- c', str(tuple(Nes_values))
    print >> out, 't1 <- c', str(tuple(t1_exact))
    print >> out, 't1.approx <- c', str(tuple(t1_approx))
    print >> out, mk_call_str('plot', 'Nes', 't1',
            #col='"red"', type='"l"', xaxp='c(0,1,10)',
            ylim='c(0,%s)' % (N_diploid*4))
    print >> out, mk_call_str('lines', 'Nes', 't1.approx', col='"red"')
    #print >> out, mk_call_str('points', 'p.20', 't1.20.exact',
            #col='"blue"')
    #
    """
    print >> out, 't1.10.approx <- c', str(tuple(t1_10_approx))
    print >> out, 't1.10.exact <- c', str(tuple(t1_10_exact))
    print >> out, 'p.10 <- c', str(tuple(proportions_10.tolist()))
    print >> out, mk_call_str('lines', 'p.10', 't1.10.approx',
            col='"red"')
    print >> out, mk_call_str('points', 'p.10', 't1.10.exact',
            col='"blue"')
    """
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

