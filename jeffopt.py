"""
This module has a maximimization function implemented by Jeff Thorne.

The maximization function has been extracted from its context by Liwen Zou.
Speed is not really important in this function.
Most of the delay should be within the evaluation itself.
"""

import unittest

import numpy as np
from scipy import special


# Stop likelihood maximization when likelihood can no longer be improved  
# by at least STOPTHRESH in a cycle.
STOPTHRESH = 1e-6
#STOPTHRESH = 1e-3

# what proportion of parameter value to be used when calculating approx. slope
INCFRAC = 0.01
#INCFRAC = 0.1

# multiply slope by this if 2nd deriv > 0  
STEPFRAC = 0.1
#STEPFRAC = 0.5

def _fmax_jeff(f, X, args, i, oldparam, incstep, lasthood):
    """
    Maximize the function value for a single parameter.
    The parameter X[i] is modified in-place.
    @param f: objective function
    @param X: current parameter values
    @param args: extra args for the objective function
    @param i: index of parameter being optimized
    @param oldparam: old value of the single parameter being optimized
    @param incstep: step size and direction
    @param lasthood: last likelihood
    @return: objective function value lhood
    """
    while True:
        if oldparam + incstep > 1 or oldparam + incstep < 0:
            incstep /= 2
        else:
            X[i] = oldparam + incstep
            lhood = f(X, *args)
            print lhood
            if lhood >= lasthood:
                return lhood
            elif X[i] == oldparam:
                # If no improvement, stop.
                return lhood
            else:
                # Newton step reduced by half
                # in order to find a point (x, x+newton-step)
                # which has largest lk.
                # New point x + newton-step/2^m m is the number of iteration
                # in this while loop.
                incstep /= 2

def fmax_jeff(f, X_guess, args=()):
    """
    Parameters have values between 0 and 1.
    Initial guesses are provided.
    Some amount of effort has been taken to match
    the minimization signatures in scipy.
    @param f: objective function
    @param X_guess: initial guess of parameters
    @param args: extra args for the objective function
    @return: best guess, best objective function value
	"""
    X = X_guess.copy()
    nparams = len(X)
    #
    oldhood = None
    lhood = None
    #
    while (oldhood is None or lhood - oldhood > STOPTHRESH):
        oldhood = lhood
        for i in range(nparams):
            lasthood = lhood
            oldparam = X[i]
            inc = INCFRAC * min((X[i], 1-X[i]))
            while oldparam + inc > 1 or oldparam - inc < 0:
                inc /= 2
            if lasthood is None:
                nowhood = f(X, *args)
                lasthood = nowhood
            else:
                nowhood = lasthood
            if not inc:
                return nowhood, X
            X[i] += inc
            plushood = f(X, *args)
            X[i] -= 2*inc
            minushood = f(X, *args)
            if plushood + minushood - 2*nowhood < 0:
                # newton step
                incstep = -0.5 * inc * (plushood - minushood) / (
                        plushood + minushood - 2*nowhood)
            else:
                # little slope step
                incstep = STEPFRAC * (plushood - nowhood) / inc
            lhood = _fmax_jeff(f, X, args, i, oldparam, incstep, lasthood)
    return lhood, X

class NegWrap:
    def __init__(self, f):
        self.f = f
    def __call__(self, *args):
        return -self.f(*args)

class LogitWrap:
    def __init__(self, f):
        self.f = f
    def __call__(self, X, *args):
        return self.f(special.logit(X), *args)

def fmin_jeff(f, X_guess, args=()):
    return fmax_jeff(NegWrap(f), X_guess, args)

def fmin_jeff_unconstrained(f, X_guess, args=()):
    lhood, X = fmin_jeff(LogitWrap(f), special.expit(X_guess), args)
    return lhood, special.logit(X)

def _ftest(X):
    # This is an upward curving parabola with min at 0.5.
    x = X[0]
    return x**2 - x + 0.25

class TestOpt(unittest.TestCase):
    def test_fmax_jeff(self):
        X_guess = np.array([0.123])
        print 'fmax'
        print X_guess
        print fmax_jeff(_ftest, X_guess)
    def test_fmin_jeff(self):
        X_guess = np.array([0.123])
        print 'fmin'
        print X_guess
        print fmin_jeff(_ftest, X_guess)

