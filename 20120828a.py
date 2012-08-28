"""
Reproduce a figure from a paper by Christina Chen et al.

"Effects of dominance on the probability of fixation of a mutant allele."
Christina T. L. Chen and Quo-Shin Chi and Stanley A. Sawyer.
2008, J. Math. Biol.
Compute fixation probability as a function
of population size, selection strength, dominance effect,
and initial allele frequency.
Fitnesses are
AA: 1 + sigma
aA: 1 + h * sigma
aa: 1
The allele 'a' is the background allele.
The scaled selection coefficient is
s = 2*N*sigma
where N is the effective population size.
To compute the fixation probability we first transform the parameters
and then we break the computation into two cases.
alpha = h / (2h - 1)
beta = 2h - 1
With this transformation, beta*s > 0 corresponds to overdominant selection
and beta*s < 0 corresponds to underdominant selection.
Let L = sqrt(abs(beta*s))
then we can compute the fixation probability as the ratio of two
definite integrals of somewhat simple functions
that depend on the sign of beta*s.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import special

import Form
import FormOut
import RUtil
from RUtil import mk_call_str

def erfi(x):
    """
    Why does scipy.special not include this?
    """
    return 2 * x * special.hyp1f1(0.5, 1.5, x*x) / math.sqrt(math.pi)

def get_form():
    """
    @return: the body of a form
    """
    return [
            #Form.Integer('N_diploid', 'diploid population size',
                #5, low=3, high=40),
            Form.ImageFormat(),
            ]

def get_fixation_probability(p, s, h):
    """
    Assume nontrivial deviation from neutrality and selective additivity.
    @param p: initial allele frequency
    @param s: positive when the mutant allele is fitter
    @param h: dominance
    @return: fixation probability
    """
    # TODO add a special case for selective additivity.
    beta = 2.0 * h - 1.0
    if not s:
        return p
    if not beta:
        raise ValueError('additive selection is not implemented')
    alpha = h / beta
    if beta * s > 0:
        # overdominant if 0 < alpha < 1
        f = erfi
    elif beta * s < 0:
        # underdominant if 0 < alpha < 1
        f = special.erf
    else:
        raise ValueError
    L = math.sqrt(abs(beta * s))
    a0 = f(L*(0 - alpha))
    a1 = f(L*(p - alpha))
    b0 = f(L*(0 - alpha))
    b1 = f(L*(1 - alpha))
    pfix = (a1 - a0) / (b1 - b0)
    if not pfix:
        print 'p:', p
        print 's:', s
        print 'h:', h
        print 'beta:', beta
        print 'alpha:', alpha
        print 'beta*s:', beta*s
        print 'pfix:', pfix
        print 'L*(0 - alpha)', L*(0 - alpha)
        print 'L*(p - alpha)', L*(p - alpha)
        print 'L*(0 - alpha)', L*(0 - alpha)
        print 'L*(1 - alpha)', L*(1 - alpha)
        print 'a0:', a0
        print 'a1:', a1
        print 'b0:', b0
        print 'b1:', b1
        print
    return pfix

def get_form_out():
    return FormOut.Image('plot')

def get_response_content(fs):
    # define the initial frequency
    p = 0.4
    # define some selection coefficients to plot
    #low = -2
    low = -10
    #high = 5
    high = 100
    s_values = np.linspace(low, high, (high-low)*10 + 1)
    # define some dominance coefficients
    h_values = [-0.5, -2.0, -3.0]
    # get the values for each h
    arr = []
    for h in h_values:
        arr.append([get_fixation_probability(p, s, h) for s in s_values])
    #
    # define the r script
    out = StringIO()
    print >> out, 'title.string <- "my title"'
    print >> out, 's.values <- c', str(tuple(s_values))
    print >> out, 'ha <- c', str(tuple(arr[0]))
    print >> out, 'hb <- c', str(tuple(arr[1]))
    print >> out, 'hc <- c', str(tuple(arr[2]))
    print >> out, mk_call_str('plot', 's.values', 'ha',
            type='"l"',
            xlab='"selection coefficient (s)"',
            ylab='"probability of ultimate fixation"',
            main='"fixation probabilities for various h ; p0 = 0.4"',
            ylim='c(0,1)',
            )
    print >> out, mk_call_str('lines', 's.values', 'hb', col='"red"')
    print >> out, mk_call_str('lines', 's.values', 'hc', col='"green"')
    script = out.getvalue().rstrip()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter_no_table(
            script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

