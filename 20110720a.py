"""
Examine properties of the Chebyshev polynomials. [UNFINISHED]

In particular look at their orthogonality,
their interlacing properties, their trigonometric relations,
and their connection to eigenfunctions of Sturm-Liouville systems.
"""

from StringIO import StringIO

import numpy as np
import scipy
from scipy import linalg
import sympy

import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    out = StringIO()
    polys = [sympy.polys.orthopolys.chebyshevt_poly(
        i, sympy.abc.x) for i in range(4)]
    for i, p in enumerate(polys):
        if i:
            print >> out, i
            print >> out, p
            print >> out, p(sympy.cos(sympy.abc.x))
            print >> out
    return out.getvalue()

