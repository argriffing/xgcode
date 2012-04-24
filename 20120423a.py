"""
Try to reproduce an attempt to calculate mutual information.

This is a variant of mutual information that is umm creative.
"""

from StringIO import StringIO
import argparse
import math

import numpy as np

import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Float('t', 'time', '0.5', low_exclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.Report('summary')

def get_conditional_selection(N, M, t):
    a = (M - 1) / M
    b = math.exp(-t*M/(N-1))
    c = 1 / M
    return a*b + c

def get_conditional_neutral(N, t):
    return get_conditional_selection(N, N, t)

def get_expectation_inf_neutral(N):
    return 1 / (N*N)

def get_expectation_inf_selection(N, M):
    return 1 / (M*M)

def get_expectation_t_neutral(N, t):
    p = get_conditional_neutral(N, t)
    a = (1 / (N-1)) * p * p
    b = (2 / (N*(N-1))) * p
    c = 1 / (N*(N-1))
    return a - b + c

def get_expectation_t_selection(N, M, t):
    p = get_conditional_selection(N, M, t)
    a = (1 / (M-1)) * p * p
    b = (2 / (M*(M-1))) * p
    c = 1 / (M*(M-1))
    return a - b + c

def get_expectation_t_neutral_tweaked(N, t):
    p = get_conditional_neutral(N, t)
    a = p * p
    b = (2 * (N-1) * p) / N
    c = (N - 1) / N
    return a - b + c

def get_expectation_t_selection_tweaked(N, M, t):
    p = get_conditional_selection(N, M, t)
    a = p * p
    b = (2 * (M-1) * p) / M
    c = (M - 1) / M
    return a - b + c


def get_response_content(fs):
    N = 4.0
    M = 2.0
    t = fs.t
    result_neutral = sum([
        math.log(get_expectation_t_neutral(N, t)),
        -math.log(get_expectation_inf_neutral(N))])
    result_selection = sum([
        math.log(get_expectation_t_selection(N, M, t)),
        -math.log(get_expectation_inf_selection(N, M))])
    result_neutral_tweaked = sum([
        math.log(get_expectation_t_neutral_tweaked(N, t)),
        -math.log(get_expectation_inf_neutral(N))])
    result_selection_tweaked = sum([
        math.log(get_expectation_t_selection_tweaked(N, M, t)),
        -math.log(get_expectation_inf_selection(N, M))])
    out = StringIO()
    print >> out, 't:', t
    print >> out
    print >> out, 'observed results (wrong):'
    print >> out, 'N=4:', result_neutral
    print >> out, 'N=4, M=2:', result_selection
    print >> out
    print >> out, 'tweaked results (extra wrong):'
    print >> out, 'N=4:', result_neutral_tweaked
    print >> out, 'N=4, M=2:', result_selection_tweaked
    return out.getvalue()

