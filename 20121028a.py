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
import matplotlib.pyplot as plt
import matplotlib

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
    if not numpy.isfinite(z):
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

def do_integration_demo():
    # http://matplotlib.org/examples/pylab_examples/layer_images.html
    dc, dd = 0.001, 0.0001
    c = np.arange(1e-16, 0.1, dc)
    d = np.arange(-0.05, 0.05, dd)
    C, D = np.meshgrid(c, d)
    cmin, cmax, dmin, dmax = np.amin(c), np.amax(c), np.amin(d), np.amax(d)
    extent = cmin, cmax, dmin, dmax
    fig = plt.figure()
    #Z = np.vectorize(kimrecessive.denom_quad)(C, D)

    #print C
    #print D

    Z = np.vectorize(get_relative_error_a)(C, D)
    im = plt.imshow(
            #Z, cmap=plt.cm.jet, interpolation='bilinear', extent=extent)
            Z, cmap=plt.cm.jet, extent=extent)
    plt.show()

    Z = np.vectorize(get_relative_error_b)(C, D)
    im = plt.imshow(
            #Z, cmap=plt.cm.jet, interpolation='bilinear', extent=extent)
            Z, cmap=plt.cm.jet, extent=extent)
    plt.show()

    """
    for c in numpy.linspace(-3, 3, 11):
        for d in numpy.linspace(-0.05, 0.05, 21):
            x = denom_piecewise(c, d)
            y = denom_quad(c, d)
            z = d**2 + (d/c)**2
            print 'c:         ', c
            print 'd:         ', d
            print 'quad:      ', y
            print 'piecewise: ', x
            print 'method:    ', z
            print denom_not_genic(c, d)
            print denom_near_genic(c, d)
            if abs(y - x) / y < 1e-6:
                print 'ok'
            else:
                print '*** bad ***'
            print
    raise Exception
    """

def main():
    do_integration_demo()

if __name__ == '__main__':
    main()

