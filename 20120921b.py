"""
Check confounding of h and s using the notation of Christina Chen et al.

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
"""

from StringIO import StringIO
import math
import cgi

import numpy as np
from scipy import linalg

import Form
import FormOut
import RUtil
from RUtil import mk_call_str
import MatrixUtil
import StatsUtil
import kimura
import wfengine

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size',
                100, low=3, high=200),
            Form.Integer('nmutants', 'initial number of mutant alleles',
                80, low=1, high=400),
            Form.Sequence('s_values', 'selection values',
                ('10', '20', '30')),
            Form.Sequence('h_values', 'dominance values',
                ('-0.5', '-2', '-3')),
            ]

def get_presets():
    return [
            Form.Preset(
                'fig 1',
                {
                    'N_diploid' : '100',
                    'nmutants' : '20',
                    's_values' : ('0', '2.5', '5'),
                    'h_values' : ('-3', '0', '1', '2'),
                    }),
            Form.Preset(
                'fig 2',
                {
                    'N_diploid' : '100',
                    'nmutants' : '160',
                    's_values' : ('-1', '0', '1'),
                    'h_values' : ('2', '3', '5', '7'),
                    }),
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

def get_pfix_finite(P):
    """
    Use the theory of absorbing Markov chains.
    The first state of P means fixation of allele 1.
    The last state of P means fixation of allele 0.
    @param P: a transition matrix
    @return: vector of fixation probabilities
    """
    n = len(P)
    Q = P[1:-1, 1:-1]
    R = P[1:-1, -1]
    I = np.eye(n-2)
    B = linalg.solve(I - Q, R)
    return B

def get_response_content(fs):
    out = StringIO()
    #
    N_hap = 2 * fs.N_diploid
    if fs.nmutants >= N_hap:
        raise Exception('too many initial mutant alleles')
    p0 = fs.nmutants / float(N_hap)
    h_values = [float(h) for h in fs.h_values]
    s_values = [float(s) for s in fs.s_values]
    #
    print >> out, '<html><body><table border="1" cellpadding="10">'
    #
    headers = (
            'method', 'h', 's',
            'fAA', 'faA', 'faa',
            'pfix',
            )
    print >> out, get_html_table_row(headers)
    #
    for h in h_values:
        for s in s_values:
            for method_name in ('finite', 'diffusion'):
                sigma = s / float(N_hap)
                #
                fAA = 1.0 + sigma
                faA = 1.0 + h * sigma
                faa = 1.0
                #
                if method_name == 'diffusion':
                    pfix = kimura.get_fixation_probability_chen(p0, s, h)
                elif method_name == 'finite':
                    P = np.exp(wfengine.create_diallelic_chen(
                        fs.N_diploid, fAA, faA, faa))
                    MatrixUtil.assert_transition_matrix(P)
                    v = get_pfix_finite(P)
                    pfix = v[fs.nmutants]
                else:
                    raise Exception('internal error')
                values = (
                        method_name, h, s,
                        fAA, faA, faa,
                        pfix,
                        )
                print >> out, get_html_table_row(values)
    print >> out, '</table></body></html>'
    return out.getvalue()

