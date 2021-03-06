"""
Check subdomains for numerical stability of a Kimura integral solution. [UNFINISHED]

Right now this only works from the command line not the internet.
"""

from StringIO import StringIO
import math
import cmath
import random

import numpy
import numpy as np
import scipy
import scipy.special
import scipy.integrate
import algopy

import kimengine
import kimrecessive
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')


def get_response_content(fs):
    #
    pass

def refilter(z):
    """
    @param z: a relative error
    @return: a more qualitative description of the error
    """
    if numpy.isnan(z):
        return 7.
    elif not numpy.isfinite(z):
        return 6.
    elif abs(z) < 1e-12:
        return 0.
    elif abs(z) < 1e-8:
        return 1.
    elif abs(z) < 1e-4:
        return 2.
    elif abs(z) < 1e0:
        return 3.
    elif abs(z) < 1e4:
        return 4.
    else:
        return 5.

def get_relative_error_a(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_not_genic(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_b(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_near_genic(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_c(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_piecewise(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_d(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_near_genic_combo(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_e1(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_hyp1f1_b(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_e2(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_hyperu_b(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_e3(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_hyp2f0_b(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_e4(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_dawsn_b(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_e5(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_erfi_b(c, d)
    z = (y-x) / x
    return refilter(z)

def get_relative_error_e6(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_combo_b(c, d)
    z = (y - x) / x
    return refilter(z)

def get_relative_error_e7(c, d):
    x = kimrecessive.denom_quad(c, d)
    y = kimrecessive.denom_poly_b(c, d)
    z = (y - x) / x
    return refilter(z)

def do_integration_demo():
    N = 101
    #d = np.linspace(-5, 5, N) / 10000.
    #c = np.linspace(-100, 100, N) / 10.
    #d = np.linspace(-5, 5, N) / 3.
    #c = np.linspace(-100, 100, N) / 3.
    d = np.linspace(-5, 5, N)
    c = np.linspace(-100, 100, N)
    #d = np.linspace(-1, 1, N) * 0.25
    #c = np.linspace(-1, 1, N) * 0.02
    #d = np.linspace(-3, 3, N) / 10.
    #c = np.linspace(-30, 30, N) / 10.
    #d = np.linspace(-3, 3, N) / 10000.
    #c = np.linspace(-30, 30, N) / 10000.
    #d = np.linspace(-3, 3, N)
    #c = np.linspace(-30, 30, N)
    #d = np.linspace(-4, 4, N)
    #c = np.linspace(-40, 40, N)
    #c = np.linspace(-20, 20, N)
    #d = np.linspace(-.2, .2, N)
    #c = np.linspace(-100, 100, N)
    #c = np.linspace(-200, 200, N)
    #dc, dd = 0.001, 0.001
    #c = np.arange(1.-0.05, 1.+0.05, dc)
    #d = np.arange(1.-0.05, 1.+0.05, dd)
    ##dc, dd = 0.005, 0.005
    ##c = np.arange(0, 0.3, dc)
    ##d = np.arange(-0.2, 0.2, dd)


    """
    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_a(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()
    """

    """
    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_b(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()
    """

    """
    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_c(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()
    """

    """
    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_d(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()
    """

    """
    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_e1(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()

    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_e2(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()

    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_e3(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()

    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_e4(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()

    Z = np.zeros((len(d), len(c)))
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            Z[j, i] = get_relative_error_e6(ci, dj)
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()
    """

    Z = np.zeros((len(d), len(c)))
    W = np.zeros((len(d), len(c)))
    quad_x, quad_w = kimrecessive.precompute_quadrature(0, 1, 101)
    for j, dj in enumerate(d):
        for i, ci in enumerate(c):
            x = kimrecessive.denom_quad(ci, dj)
            #y = kimrecessive.denom_poly_b(ci, dj)
            #y = kimengine.denom_poly(ci, dj)
            #y = kimrecessive.denom_hyperu_b(ci, dj)
            #y = kimrecessive.denom_erfcx_b(ci, dj)
            y = kimrecessive.denom_fixed_quad(ci, dj, quad_x, quad_w)
            #y = kimrecessive.denom_combo_b(ci, dj)
            w = abs(y - x) / abs(x)
            W[j, i] = w
            Z[j, i] = refilter(w)
    print numpy.max(W)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    im = plt.imshow(Z, cmap=plt.cm.jet)
    plt.show()


def main():
    do_integration_demo()

if __name__ == '__main__':
    main()

