"""
Check the generalized fixation probability using the notation of Kimura.

The generalization is to allow recessiveness/dominance;
the previous notation was that of Chen et al. which was less standard.
When the dominance h value is 0.5,
the formulas in this web script should correspond exactly to Kimura's
genic approximations of fixation probability
of a single mutant allele
with selection and drift but no mutational pressure,
recombination, or migration.
In this web script,
the "effective population size" is assumed to be equal to the
diploid population size.
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
            Form.Integer('N', 'Kimura diploid population N',
                10000, low=1),
            Form.Sequence('x4Ns_values', 'Kimura 4Ns selection values',
                ('1', '2', '4')),
            Form.Sequence('h_values', 'Kimura(?) dominance h values',
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

def get_pfix_transformed(p0, x4Ns, h):
    """
    Try to get the same result as the function in the kimura module.
    Change of variables according to eqn 3 of Chen et al.
    When I type
    (integral from 0 to p of exp(b*s*(x-a)**2) ) /
    (integral from 0 to 1 of exp(b*s*(x-a)**2) )
    I get
    ( erfi(sqrt(b)*sqrt(s)*a) - erfi(sqrt(b)*sqrt(s)*(a-p)) ) /
    ( erfi(sqrt(b)*sqrt(s)*a) - erfi(sqrt(b)*sqrt(s)*(a-1)) )
    @param p0: proportion of mutant alleles in the population
    @param x4Ns: 4Ns
    @return: fixation probability
    """
    if not x4Ns:
        # This is the neutral case.
        return p0
    if h == 0.5:
        # This is the genic case.
        # Checking for exact equality of 0.5 is OK.
        return math.expm1(-x4Ns*p0) / math.expm1(-x4Ns)
    b = 2.0 * h - 1.0
    a = h / (2.0 * h - 1.0)
    q = cmath.sqrt(b) * cmath.sqrt(x4Ns)
    #
    top = kimura.erfi(q*a) - kimura.erfi(q*(a-p0))
    bot = kimura.erfi(q*a) - kimura.erfi(q*(a-1))
    return top / bot

def get_pfix_transformed_limit(x4Ns, h):
    """
    This is an approximation analogous to an approximation by Kimura.
    @param x4Ns: 4Ns
    @param h: dominance parameter
    @return: 2N times fixation probability
    """
    if not x4Ns:
        # This is the neutral case.
        return 1.0
    if h == 0.5:
        # This is the genic case.
        # Checking for exact equality of 0.5 is OK.
        return - x4Ns / math.expm1(-x4Ns)
    a = h / (2.0 * h - 1.0)
    b = 2.0 * h - 1.0
    q = cmath.sqrt(b) * cmath.sqrt(x4Ns)
    #
    top = 2*q*cmath.exp((q*a)**2)
    bot_a = kimura.erfi(q*a)
    bot_b = kimura.erfi(q*(1-a))
    return top / ( cmath.sqrt(math.pi) * (bot_a + bot_b) )

#XXX dont know what is up with this function
def get_pfix(p0, s, h):
    a = -s
    b = 2*h
    q = cmath.sqrt(a) / (2 * cmath.sqrt(b))
    top_a = special.erf(q*(b*(2*p0 - 1) - 1))
    top_b = special.erf(q*(b+1))
    bot_a = special.erf(q*(b+1))
    bot_b = special.erf(q*(b-1))
    return (top_a + top_b) / (bot_a + bot_b)

#XXX This function is wrong.
def scaled_limit(s_in, h):
    """
    This was just wrong because I used the wrong equation.
    The next line should have ...(x + b*... instead.
    ( integral from 0 to p of exp( a*x*(1 + b*(1-x)) ) ) ) /
    ( integral from 0 to 1 of exp( a*x*(1 + b*(1-x)) ) ) )
    The first term of the series expansion of the integral
    around p=0, divided by p.
    Note that p can be something like 1/(2*N).
    Dividing by p in this first term of the integral expansion
    removes the dependence on p.
    """
    s = 2*s_in
    a = -s
    b = 2*h
    q = cmath.sqrt(a) / (2 * cmath.sqrt(b))
    top = 2*cmath.sqrt(a)*cmath.sqrt(b)*cmath.exp(-(b+1)*(b+1)*a/(4*b))
    bot_a = special.erf(q*(b+1))
    bot_b = special.erf(q*(b-1))
    return top / ( cmath.sqrt(math.pi) * (bot_a + bot_b) )

#XXX This is probably also wrong.
def second_order_approx(p0, s, h):
    a = -s
    b = 2*h
    term_1 = p0 * scaled_limit(s, h)
    term_2 = 0.5 * a * (b+1) * p0 * term_1
    return term_1 + term_2

def get_response_content(fs):
    out = StringIO()
    #
    N = fs.N
    nmutants = 1
    p0 = nmutants / float(2*N)
    h_values = [float(h) for h in fs.h_values]
    x4Ns_values = [float(x4Ns) for x4Ns in fs.x4Ns_values]
    #
    print >> out, '<html><body><table border="1" cellpadding="10">'
    #
    headers = (
            'kimura h', 'kimura 4Ns', 'kimura s',
            'w11', 'w12', 'w22',
            'diffusion exact',
            'diffusion approx',
            )
    print >> out, get_html_table_row(headers)
    #
    for h in h_values:
        for x4Ns in x4Ns_values:
            s = x4Ns / float(4*N)
            #
            w11 = 1.0 + s
            w12 = 1.0 + h * s
            w22 = 1.0
            #
            diffusion_exact = get_pfix_transformed(p0, x4Ns, h)
            diffusion_approx = p0 * get_pfix_transformed_limit(x4Ns, h)
            values = (
                    h, x4Ns, s,
                    w11, w12, w22,
                    diffusion_exact.real,
                    diffusion_approx.real,
                    )
            print >> out, get_html_table_row(values)
    print >> out, '</table></body></html>'
    return out.getvalue()

