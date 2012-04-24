"""
Check a formula for Jukes-Cantor plus selection.
"""

from StringIO import StringIO
import argparse
import math
from math import exp, log, cosh, sinh, tanh

import numpy as np
import scipy
from scipy import linalg

import Form
import FormOut
from MatrixUtil import ndot
import mrate
import ctmcmi
import RUtil


def coth(x):
    return cosh(x) / sinh(x)

def sech(x):
    return 1 / cosh(x)

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_jc_rate_matrix():
    n = 4
    R = np.ones((4, 4), dtype=float)
    R -= np.diag(np.sum(R, axis=1))
    return R

def get_mut_sel_rate_matrix(y, z):
    """
    This is a mutation-selection nucleotide rate matrix.
    It looks like
    * 1 y y
    1 * y y
    z z * 1
    z z 1 *
    """
    n = 4
    R = np.zeros((4, 4), dtype=float)
    for i in range(n):
        for j in range(n):
            if i in (0, 1) and j in (2, 3):
                R[i, j] = y
            elif i in (2, 3) and j in (0, 1):
                R[i, j] = z
            elif i != j:
                R[i, j] = 1
    R -= np.diag(np.sum(R, axis=1))
    return R

def x_to_halpern_bruno_yz(x):
    """
    @param x: ratio of probability of second group to first group
    @return: y, z
    """
    if x == 1:
        y = 1.0
        z = 1.0
    else:
        tau_ab = x
        tau_ba = 1/x
        y = math.log(tau_ab) / (1 - 1 / tau_ab)
        z = math.log(tau_ba) / (1 - 1 / tau_ba)
    return y, z

def get_trans_mat_expm(R, t):
    P = scipy.linalg.expm(R*t)
    return P

def get_trans_mat_tediously(y, z, t):
    """
    This is from MatrixExp in wolfram alpha.
    """
    e1m = math.exp(-2*t*(y+z))
    e1p = math.exp(2*t*(y+z))
    e2m = math.exp(-2*t*(y+1))
    e2p = math.exp(2*t*(y+1))
    e3m = math.exp(-2*t*(z+1))
    e3p = math.exp(2*t*(z+1))
    e4 = math.exp(2*t*(y+1+y+z))
    e5 = math.exp(2*t*(z+1+y+z))
    pre = e1m / (y+z)
    a = pre * 0.5 * e2m * (y*e2p - (y+z)*e1p + z*e4)
    b = pre * 0.5 * e3m * (z*e3p - (y+z)*e1p + y*e5)
    c = pre * 0.5 * y * (e1p - 1)
    d = pre * 0.5 * z * (e1p - 1)
    n = 4
    P = np.array([
        [0, a, c, c],
        [a, 0, c, c],
        [d, d, 0, b],
        [d, d, b, 0]])
    for i in range(n):
        P[i,i] = 1.0 - np.sum(P[i])
    return P

def get_trans_mat_tediously_b(y, z, t):
    """
    Algebraic simplification.
    """
    py = y / (y + z)
    pz = z / (y + z)
    eyz = math.exp(-2*t*(y+z))
    ey1 = math.exp(-2*t*(y+1))
    ez1 = math.exp(-2*t*(z+1))
    a = 0.5 * (py*eyz + pz - ey1)
    b = 0.5 * (pz*eyz + py - ez1)
    c = 0.5 * py * (1 - eyz)
    d = 0.5 * pz * (1 - eyz)
    n = 4
    P = np.array([
        [0, a, c, c],
        [a, 0, c, c],
        [d, d, 0, b],
        [d, d, b, 0]])
    for i in range(n):
        P[i,i] = 1.0 - np.sum(P[i])
    return P

def get_trans_mat_tediously_c(y, z, t):
    """
    Algebraic simplification.
    This time py and pz are actually stationary probabilities.
    """
    py = 0.5 * z / (y + z)
    pz = 0.5 * y / (y + z)
    eyz = math.exp(-2*t*(y+z))
    ey1 = math.exp(-2*t*(y+1))
    ez1 = math.exp(-2*t*(z+1))
    a = pz*eyz + py - 0.5 * ey1
    b = py*eyz + pz - 0.5 * ez1
    c = pz * (1 - eyz)
    d = py * (1 - eyz)
    n = 4
    P = np.array([
        [0, a, c, c],
        [a, 0, c, c],
        [d, d, 0, b],
        [d, d, b, 0]])
    for i in range(n):
        P[i,i] = 1.0 - np.sum(P[i])
    return P

def get_trans_mat_from_x(x, t):
    """
    Algebraic simplification.
    Assume Halpern-Bruno.
    Do not go through y and z.
    """
    py = 0.5 * 1 / (x + 1)
    pz = 0.5 * x / (x + 1)
    eyz = x ** (-2 * t * (x+1) / (x-1))
    ey1 = math.exp(-2*t) * ( x ** (-2*t*x / (x-1) ) )
    ez1 = math.exp(-2*t) * ( x ** (-2*t*1 / (x-1) ) )
    a = pz*eyz + py - 0.5 * ey1
    b = py*eyz + pz - 0.5 * ez1
    c = pz * (1 - eyz)
    d = py * (1 - eyz)
    n = 4
    P = np.array([
        [0, a, c, c],
        [a, 0, c, c],
        [d, d, 0, b],
        [d, d, b, 0]])
    for i in range(n):
        P[i,i] = 1.0 - np.sum(P[i])
    return P

def get_jc_e_ll(t):
    """
    An algebraic expression.
    """
    pa = 0.25 + 0.75 * exp(-4*t)
    pb = 0.25 - 0.25 * exp(-4*t)
    return pa*log(4*pa) + 3*pb*log(4*pb)

def get_e_ll_from_x(x, t):
    """
    Algebraic simplification.
    Assume Halpern-Bruno.
    Do not go through y and z.
    """
    py = 0.5 * 1 / (x + 1)
    pz = 0.5 * x / (x + 1)
    eyz = x ** (-2 * t * (x+1) / (x-1))
    ey1 = math.exp(-2*t) * ( x ** (-2*t*x / (x-1) ) )
    ez1 = math.exp(-2*t) * ( x ** (-2*t*1 / (x-1) ) )
    a = pz*eyz + py - 0.5 * ey1
    b = py*eyz + pz - 0.5 * ez1
    c = pz * (1 - eyz)
    d = py * (1 - eyz)
    a_diag = 1 - (a + 2*c)
    b_diag = 1 - (b + 2*d)
    # get some things divided by stationary probabilities
    a_dpy = eyz*x + 1 - (x+1)*ey1
    b_dpz = eyz/x + 1 - (x+1)*ez1/x
    #a_dpy = a / py
    #b_dpz = b / pz
    c_dpz = 1 - eyz
    d_dpy = 1 - eyz
    a_diag_dpy = a_diag / py
    b_diag_dpz = b_diag / pz
    # get the two KL divergences for the two distinct row types
    KL_y = sum([
        a_diag * math.log(a_diag_dpy),
        a * math.log(a_dpy),
        2 * c * math.log(c_dpz)])
    KL_z = sum([
        b_diag * math.log(b_diag_dpz),
        b * math.log(b_dpz),
        2 * d * math.log(d_dpy)])
    e_ll_ratio = 2 * py * KL_y + 2 * pz * KL_z
    return e_ll_ratio

def get_e_ll_from_x_b(x, t):
    """
    Further algebraic simplification.
    """
    py = 0.5 * 1 / (x + 1)
    pz = 0.5 * x / (x + 1)
    eyz = x ** (-2 * t * (x+1) / (x-1))
    ey1 = math.exp(-2*t) * ( x ** (-2*t*x / (x-1) ) )
    ez1 = math.exp(-2*t) * ( x ** (-2*t*1 / (x-1) ) )
    a = pz*eyz + py - 0.5 * ey1
    b = py*eyz + pz - 0.5 * ez1
    c = pz * (1 - eyz)
    d = py * (1 - eyz)
    a_diag = 1 - (a + 2*c)
    b_diag = 1 - (b + 2*d)
    # get some things divided by stationary probabilities
    a_dpy = eyz*x + 1 - (x+1)*ey1
    b_dpz = eyz/x + 1 - (x+1)*ez1/x
    #a_dpy = a / py
    #b_dpz = b / pz
    c_dpz = 1 - eyz
    d_dpy = 1 - eyz
    a_diag_dpy = a_diag / py
    b_diag_dpz = b_diag / pz
    # combine some stuff
    pypz = 0.25 * x / (x + 1)**2
    pya = pypz*eyz + py*py - 0.5 * py * ey1
    pzb = pypz*eyz + pz*pz - 0.5 * pz * ez1
    # get the two KL divergences for the two distinct row types
    a_diag_contrib = 2 * py * a_diag * math.log(a_diag_dpy)
    a_contrib = 2 * pya * math.log(a_dpy)
    b_diag_contrib = 2 * pz * b_diag * math.log(b_diag_dpz)
    b_contrib = 2 * pzb * math.log(b_dpz)
    dc_contrib = 8 * pypz * (1 - eyz) * math.log(1 - eyz)
    contrib_arr = [
        a_diag_contrib,
        a_contrib,
        b_diag_contrib,
        b_contrib,
        dc_contrib]
    print contrib_arr
    e_ll_ratio = sum(contrib_arr)
    return e_ll_ratio

def get_e_ll_from_x_htrig(x, t):
    """
    Further algebraic simplification.
    Try to use hyperbolic trigonometric functions.
    """
    w = 0.5 * math.log(x)
    py = (1 - tanh(w))/4
    pz = (1 + tanh(w))/4
    eyz = exp(-4*t*w*coth(w))
    ey1 = exp(-2*t*(w*coth(w) + w + 1))
    ez1 = exp(-2*t*(w*coth(w) - w + 1))
    #a = pz*eyz + py - 0.5 * ey1
    a = (1 + tanh(w))*exp(-4*t*w*coth(w))/4 + (1-tanh(w))/4 - exp(-2*t*(w*coth(w)+w+1))/2
    #b = py*eyz + pz - 0.5 * ez1
    b = (1 - tanh(w))*exp(-4*t*w*coth(w))/4 + (1+tanh(w))/4 - exp(-2*t*(w*coth(w)-w+1))/2
    c = pz * (1 - eyz)
    d = py * (1 - eyz)
    a_diag = 1 - (a + 2*c)
    b_diag = 1 - (b + 2*d)
    # get some things divided by stationary probabilities
    #a_dpy = eyz*x + 1 - (x+1)*ey1
    a_dpy = exp(2*w-4*t*w*coth(w)) + 1 - (exp(2*w)+1)*exp(-2*t*(w*coth(w)+w+1))
    #b_dpz = eyz/x + 1 - (x+1)*ez1/x
    b_dpz = exp(-2*w-4*t*w*coth(w)) + 1 - (exp(2*w)+1)*exp(-2*w-2*t*(w*coth(w)-w+1))
    #c_dpz = 1 - eyz
    c_dpz = 1 - exp(-4*t*w*coth(w))
    #d_dpy = 1 - eyz
    d_dpy = 1 - exp(-4*t*w*coth(w))
    a_diag_dpy = a_diag / py
    b_diag_dpz = b_diag / pz
    # combine some stuff
    pypz = (sech(w)/4)**2
    pya = pypz*eyz + py*py - 0.5 * py * ey1
    pzb = pypz*eyz + pz*pz - 0.5 * pz * ez1
    # get the two KL divergences for the two distinct row types
    a_diag_contrib = 2 * py * a_diag * math.log(a_diag_dpy)
    a_contrib = 2 * pya * math.log(a_dpy)
    b_diag_contrib = 2 * pz * b_diag * math.log(b_diag_dpz)
    b_contrib = 2 * pzb * math.log(b_dpz)
    #dc_contrib = 8 * pypz * (1 - eyz) * math.log(1 - eyz)
    dc_contrib = (sech(w)**2)*(1-exp(-4*t*w*coth(w)))*log(1-exp(-4*t*w*coth(w)))/2
    contrib_arr = [
        a_diag_contrib,
        a_contrib,
        b_diag_contrib,
        b_contrib,
        dc_contrib]
    print contrib_arr
    e_ll_ratio = sum(contrib_arr)
    return e_ll_ratio

def x_to_distn(x):
    """
    @param x: ratio of stationary probabilities
    """
    v = np.array([1, 1, x, x], dtype=float)
    return v / np.sum(v)

def get_response_content(fs):
    np.set_printoptions(linewidth=200)
    out = StringIO()
    R_jc = get_jc_rate_matrix()
    t = 0.1
    x = 1.6
    w = 0.5 * log(x)
    v = x_to_distn(x)
    R_hb_easy = mrate.to_gtr_halpern_bruno(R_jc, v)
    y, z, = mrate.x_to_halpern_bruno_yz(x)
    yz_ratio = y / z
    R_hb_tedious = get_mut_sel_rate_matrix(y, z)
    P_hb_easy = get_trans_mat_expm(R_hb_easy, t)
    P_hb_tedious = get_trans_mat_tediously(y, z, t)
    P_hb_tedious_c = get_trans_mat_tediously_c(y, z, t)
    P_hb_from_x = get_trans_mat_from_x(x, t)
    e_ll_jc = ctmcmi.get_expected_ll_ratio(R_jc, t)
    e_ll_jc_tedious = get_jc_e_ll(t)
    e_ll_hb = ctmcmi.get_expected_ll_ratio(R_hb_easy, t)
    e_ll_hb_from_x = get_e_ll_from_x(x, t)
    e_ll_hb_from_x_b = get_e_ll_from_x_b(x, t)
    e_ll_hb_from_x_htrig = get_e_ll_from_x_htrig(x, t)
    # print some values
    print >> out, 'Jukes-Cantor mutation matrix:'
    print >> out, R_jc
    print >> out
    print >> out, 'ratio of common to uncommon probabilities:'
    print >> out, x
    print >> out
    print >> out, '1/2 log ratio:'
    print >> out, w
    print >> out
    print >> out, 'fast rate:'
    print >> out, y
    print >> out
    print >> out, 'slow rate:'
    print >> out, z
    print >> out
    print >> out, 'reciprocal of fast rate:'
    print >> out, 1.0 / y
    print >> out
    print >> out, 'ratio of fast to slow rates (should be x):'
    print >> out, yz_ratio
    print >> out
    print >> out, 'mutation-selection rate matrix (easy):'
    print >> out, R_hb_easy
    print >> out
    print >> out, 'mutation-selection rate matrix (tedious):'
    print >> out, R_hb_tedious
    print >> out
    print >> out, 'time:'
    print >> out, t
    print >> out
    print >> out, 'mutation-selection transition matrix (easy):'
    print >> out, P_hb_easy
    print >> out
    print >> out, 'mutation-selection transition matrix (tedious):'
    print >> out, P_hb_tedious
    print >> out
    print >> out, 'mutation-selection transition matrix (tedious c):'
    print >> out, P_hb_tedious_c
    print >> out
    print >> out, 'mutation-selection transition matrix (from x):'
    print >> out, P_hb_from_x
    print >> out
    print >> out, 'expected Jukes-Cantor log likelihood ratio:'
    print >> out, e_ll_jc
    print >> out
    print >> out, 'expected Jukes-Cantor log likelihood ratio (tedious):'
    print >> out, e_ll_jc_tedious
    print >> out
    print >> out, 'expected mutation-selection log likelihood ratio:'
    print >> out, e_ll_hb
    print >> out
    print >> out, 'expected mutation-selection ll ratio from x:'
    print >> out, e_ll_hb_from_x
    print >> out
    print >> out, 'expected mutation-selection ll ratio from x (impl b):'
    print >> out, e_ll_hb_from_x_b
    print >> out
    print >> out, 'expected mutation-selection ll ratio from x (htrig):'
    print >> out, e_ll_hb_from_x_htrig
    print >> out
    # check some invariants
    if np.allclose(R_hb_easy, R_hb_tedious):
        print >> out, 'halpern-bruno rate matrices are equal as expected'
    else:
        print >> out, '*** halpern-bruno rate matrices are not equal!'
    if np.allclose(P_hb_easy, P_hb_tedious):
        print >> out, 'halpern-bruno transition matrices are equal as expected'
    else:
        print >> out, '*** halpern-bruno transition matrices are not equal!'
    if np.allclose(P_hb_easy, P_hb_tedious_c):
        print >> out, 'halpern-bruno transition matrices are equal as expected'
    else:
        print >> out, '*** halpern-bruno transition matrices are not equal!'
    if np.allclose(P_hb_easy, P_hb_from_x):
        print >> out, 'halpern-bruno transition matrices are equal as expected'
    else:
        print >> out, '*** halpern-bruno trans. mat. from x is not equal!'
    # return the results
    return out.getvalue()

