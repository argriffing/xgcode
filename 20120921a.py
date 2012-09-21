"""
Check whether dominance and selection are confounded, in table form.
"""

from StringIO import StringIO
import time
import math
import cgi

import numpy as np
from scipy import optimize
from scipy import linalg

import Form
import FormOut
import MatrixUtil
import StatsUtil
import kaizeng
import wfengine
import wrightfisher

def get_form():
    """
    @return: the body of a form
    """
    return [
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

def get_fixation_probabilities(P):
    """
    Use the theory of absorbing Markov chains.
    The first state of P means fixation of allele 1.
    The last state of P means fixation of allele 0.
    @param P: a transition matrix
    @return: P01 and P10 fixation probabilities
    """
    n = len(P)
    Q = P[1:-1, 1:-1]
    R = P[1:-1, (0, -1)]
    I = np.eye(n-2)
    B = linalg.solve(I - Q, R)
    #print B
    P01 = B[-1, 0]
    P10 = B[0, -1]
    return P01, P10

def get_response_content(fs):
    out = StringIO()
    print >> out, '<html><body><table border="1" cellpadding="10">'
    #
    headers = (
            'N', 's', '2Ns', 'h',
            'f00', 'f01', 'f11',
            '2Nsh', 'NP01', 'NP10')
    print >> out, get_html_table_row(headers)
    # In this web script, N is a diploid population size.
    # get the transition matrix
    for x2Nsh in (0.075, 0.75):
        for h in (0.25, 0.5, 0.75):
            for N in (50, 200):
                x2N = 2*N
                x2Ns = x2Nsh / float(h)
                s = x2Ns / float(x2N)
                #
                P = np.exp(wfengine.create_diallelic_recessive(N, s, 1-h))
                #print P
                MatrixUtil.assert_transition_matrix(P)
                P01, P10 = get_fixation_probabilities(P)
                NP01 = N * P01
                NP10 = N * P10
                #
                f00 = 1.0
                f11 = 1.0 - s
                f01 = (1-h) * f00 + h * f11
                #
                values = (
                        N, s, x2Ns, h,
                        f00, f01, f11,
                        x2Nsh, NP01, NP10,
                        )
                print >> out, get_html_table_row(values)
    print >> out, '</table></body></html>'
    return out.getvalue()


