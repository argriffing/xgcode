"""
Check selection coefficients estimated within a mutation-selection model. [UNFINISHED]

This currently uses hardcoded max liklihood parameter estimates
from a model with two parameters for recessivity as a function of selection.
The idea is to check the implicit selection estimates as an intermediate step
towards checking the recessivity values implied by the combination of the
recessivity parameter estimates with the estimated selection coefficients.
"""

from StringIO import StringIO
import math

import numpy
import scipy
import scipy.special

import Form
import FormOut
import kimrecessive

def get_form():
    """
    @return: the body of a form
    """
    return [
            ]

def get_form_out():
    return FormOut.Report()


def get_response_content(fs):
    numpy.set_printoptions(linewidth=200)
    out = StringIO()
    #
    #
    return out.getvalue()

