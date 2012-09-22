"""
Check confounding of h and s using diffusion theory for large populations.

The confounding is checked with respect to fixation probability.
This script assumes a large background population size
with a single mutant allele and no subsequent mutation.
It uses the notation in the paper
"Effects of dominance on the probability of fixation of a mutant allele."
Christina T. L. Chen and Quo-Shin Chi and Stanley A. Sawyer.
2008, J. Math. Biol.
Fitnesses are
AA: 1 + sigma
aA: 1 + h * sigma
aa: 1
The allele 'a' is the background allele.
The scaled selection coefficient is
s = 2*N*sigma
where N is the effective population size.
"""

from StringIO import StringIO
import cgi
import math
import cmath
from scipy import special

import Form
import FormOut
import kimura

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size',
                10000, low=1),
            Form.Sequence('s_values', 'selection values',
                ('1', '2', '4')),
            Form.Sequence('h_values', 'dominance values',
                ('0.5', '0.25')),
            ]

def get_form_out():
    return FormOut.Html()

def get_html_table_row(arr):
    out = StringIO()
    print >> out, '<tr>'
    for value in arr:
        print >> out, '<td>%s</td>' % cgi.escape(str(value))
    print >> out, '</tr>'
    return out.getvalue().rstrip()

def get_pfix_transformed(p0, s_in, h):
    """
    Try to get the same result as the function in the kimura module.
    Change of variables according to eqn 3 of Chen et al.
    When I type
    (integral from 0 to p of exp(b*s*(x-a)**2) ) /
    (integral from 0 to 1 of exp(b*s*(x-a)**2) )
    I get
    ( erfi(sqrt(b)*sqrt(s)*a) - erfi(sqrt(b)*sqrt(s)*(a-p)) ) /
    ( erfi(sqrt(b)*sqrt(s)*a) - erfi(sqrt(b)*sqrt(s)*(a-1)) )
    """
    s = s_in * 2
    b = 2.0 * h - 1.0
    #XXX
    if not b:
        return float('nan')
    # transform h to alpha and beta, represented here by a and b
    a = h / (2.0 * h - 1.0)
    # intermediate parameter
    q = cmath.sqrt(b) * cmath.sqrt(s)
    #
    """
    top = 2*q*cmath.exp(a*a*b*s)
    bot_a = kimura.erfi(q*(1-a))
    bot_b = kimura.erfi(q*a)
    """
    top = kimura.erfi(q*a) - kimura.erfi(q*(a-p0))
    bot = kimura.erfi(q*a) - kimura.erfi(q*(a-1))
    return top / bot

def get_pfix_transformed_limit(s_in, h):
    s = s_in * 2
    b = 2.0 * h - 1.0
    #XXX
    if not b:
        return float('nan')
    # transform h to alpha and beta, represented here by a and b
    a = h / (2.0 * h - 1.0)
    # intermediate parameter
    q = cmath.sqrt(b) * cmath.sqrt(s)
    #
    top = 2*q*cmath.exp(a*a*b*s)
    bot_a = kimura.erfi(q*a)
    bot_b = kimura.erfi(q*(1-a))
    return top / ( cmath.sqrt(math.pi) * (bot_a + bot_b) )

def get_pfix(p0, s, h):
    a = -s
    b = 2*h
    q = cmath.sqrt(a) / (2 * cmath.sqrt(b))
    top_a = special.erf(q*(b*(2*p0 - 1) - 1))
    top_b = special.erf(q*(b+1))
    bot_a = special.erf(q*(b+1))
    bot_b = special.erf(q*(b-1))
    return (top_a + top_b) / (bot_a + bot_b)

def scaled_limit(s, h):
    """
    ( integral from 0 to p of exp( a*x*(1 + b*(1-x)) ) ) ) /
    ( integral from 0 to 1 of exp( a*x*(1 + b*(1-x)) ) ) )
    The first term of the series expansion of the integral
    around p=0, divided by p.
    Note that p can be something like 1/(2*N).
    Dividing by p in this first term of the integral expansion
    removes the dependence on p.
    """
    a = -s
    b = 2*h
    q = cmath.sqrt(a) / (2 * cmath.sqrt(b))
    top = 2*cmath.sqrt(a)*cmath.sqrt(b)*cmath.exp(-(b+1)*(b+1)*q*q)
    bot_a = special.erf(q*(b+1))
    bot_b = special.erf(q*(b-1))
    bot = cmath.sqrt(math.pi)*(bot_a + bot_b)
    return top / bot

def second_order_approx(p0, s, h):
    a = -s
    b = 2*h
    term_1 = p0 * scaled_limit(s, h)
    term_2 = 0.5 * a * (b+1) * p0 * term_1
    return term_1 + term_2

def get_response_content(fs):
    out = StringIO()
    #
    nmutants = 1
    N_hap = 2 * fs.N_diploid
    if nmutants >= N_hap:
        raise Exception('too many initial mutant alleles')
    p0 = nmutants / float(N_hap)
    h_values = [float(h) for h in fs.h_values]
    s_values = [float(s) for s in fs.s_values]
    #
    print >> out, '<html><body><table border="1" cellpadding="10">'
    #
    headers = (
            'h', 's',
            'fAA', 'faA', 'faa',
            'pfix',
            'first order',
            'second order',
            'inf order',
            )
    print >> out, get_html_table_row(headers)
    #
    for h in h_values:
        for s in s_values:
            sigma = s / float(fs.N_diploid)
            #
            fAA = 1.0 + sigma
            faA = 1.0 + h * sigma
            faa = 1.0
            #
            pfix = kimura.get_fixation_probability_chen(p0, s, h)
            #lim = p0 * scaled_limit(s, h)
            lim = p0 * get_pfix_transformed_limit(s, h)
            values = (
                    h, s,
                    fAA, faA, faa,
                    pfix,
                    lim,
                    float('nan'),
                    #second_order_approx(p0, s, h),
                    get_pfix_transformed(p0, s, h),
                    )
            print >> out, get_html_table_row(values)
    print >> out, '</table></body></html>'
    return out.getvalue()

