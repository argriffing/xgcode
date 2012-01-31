"""
This module is related to Taylor expansion of CTMC mutual information.
"""

import unittest
import math

import numpy as np
import scipy
from scipy import optimize
import mpmath as mp
from mpmath import mpf
from mpmath import calculus
from mpmath.calculus import optimization


def get_key_time_points(lam, p, N):
    """
    After this time point, a uniform approximation can be used for H entries.
    H is a matrix of taylor approximation errors.
    The relevant eigenvalue should be negative,
    and the minimum stationary probability should be positive.
    @param lam: the relevant eigenvalue of the rate matrix
    @param p: the minimum stationary probability
    @param N: the number of states in the rate matrix
    @return: two key time points
    """
    if not lam < 0:
        raise ValueError('the eigenvalue should be negative')
    if not p > 0:
        raise ValueError('the min stationary probability should be positive')
    # We want each Xij to be small in absolute value
    # compared to the minimum stationary probability.
    # Half the size should work for our purposes.
    # We can bound Xij using the eigenvalue.
    # X = [e^S]_{i,j} - \sqrt{pi_i pi_j}
    # X_{i,j} = \sum_k^N U_{i,k} U_{j,k} e^{\lamda_k t}
    # e^{\lambda_2 t} \leq |X_{i,j}| \leq (N-1)e^{\lambda_k t}
    # After this amount of time,
    # we have a uniform error bound on the entries of H.
    time_to_uniformity = (1 / lam) * math.log(p / (2*(N-1)))
    # Define the bounding slope of the h bound
    # for times after the time to uniformity.
    # After this time point,
    # x_abs_upper_bound * M gives an upper bound
    # on the entries of the H matrix.
    ### y_low = p
    ### x_abs_max = p/2
    ### x_low = -x_abs_max
    ### M = mi_taylor_h(x_low, y_low) / x_abs_max
    ### M = get_M(p)
    # Now we need to use M to determine a second time point.
    # This second time point will guarantee a lower bound
    # on the mutual information which excludes zero.
    # As noted in my notes jan 27 2012 [b] this inequality never holds
    ### if M*(N**2)*((N-1)**3)*math.exp(lam * time_to_uniformity) < 0.25:
        ### time_to_usefulness = time_to_uniformity
    ### else:
        ### time_to_usefulness = -(1/lam)*math.log(4*M*(N**2)*((N-1)**3))
    # define the relaxation time
    tau = - 1 / lam
    c = 4*(3 - 2*math.log(4))
    time_to_usefulness =  tau * math.log((c/p)*(N**2)*((N-1)**3))
    # At every time point beyond the time to usefulness,
    # (1/4)exp(2 lambda_2 t) <= MI(t) <= (1/2)(N-(1/2))*exp(2 lamda_2 t) .
    # Therefore this provides the necessary bounds
    # because the mutual information is bounded away from zero
    # and its upper and lower bounds have the same decay exponent.
    #
    #
    return time_to_uniformity, time_to_usefulness

def get_sophisticated_time_bound(r1, r2, N1, N2, p1, p2):
    """
    This will use the O(x^3) rather than a rougher O(x^2) error.
    These error asymptotics follow from a second order Taylor approximation.
    This sophisticated time bound will only work after
    the point at which we have a uniform bound on
    the entries of a certain matrix
    so that we can control the sum of the entries.
    We are calling this prerequisite time point the time to uniformity.
    @param r1: absolute value of the slow second eigenvalue
    @param r2: absolute value of the fast second eigenvalue
    @param N1: number of states in slow process
    @param N2: number of states in fast process
    @param p1: min stationary probability in the slow process
    @param p2: min stationary probability in the fast process
    """
    # use four times as much precision as usual
    #mp.prec = 53 * 4
    #
    # Now define a new time point.
    # This should be after uniformity,
    # and it should be computed by numerical approximation.
    # It is the intersection of two functions of time,
    # where each function is the sum of two
    # not-necessarily-positive exponentials.
    # A key substitution will be x = exp(-t).
    #
    # Let 0 < r1 < r2 be the negatives of the second eigenvalues.
    # Let 1 < N1, 1 < N2 be the integer numbers of states.
    # Let 0 < p1 < 1, 0 < p2 < 1 be the min stationary probabilities.
    #
    # The dynamics of the mutual information can be divided
    # into several time scales.
    # At the smallest scale at times near zero,
    # the mutual information is defined entirely by the
    # entropy of the stationary distribution.
    # At slightly larger time scales
    # it is determined by complicated dynamics
    # involving the stationary probabilities.
    # At even larger time scales, it becomes possible
    # to obtain bounds using relatively simple functions,
    # in particular the sum of two exponential functions.
    # This is the time scale of interest here.
    # At even larger time scales which we do not consider here,
    # it is possible to obtain informative bounds
    # that are pure exponential functions rather than the
    # sum of two exponential functions.
    #
    # The simple form is:
    # a1 exp(-2 r1 t) - b1 exp(-3 r1 t) = a2 exp(-2 r2 t) + b2 exp(-3 r2 t)
    # where the lhs is a lower bound of a slowly decaying function
    # and the rhs is an upper bound of a quickly decaying function
    #
    # Substituting x = exp(-t):
    # a1 x^(2 r1) - b1 x^(3 r1) = a2 x^(2 r2) + b2 x*(3 r2)
    #
    # define coefficient bounds of the cubic error term
    b1 = mpf( (N1**2) * ((N1-1)**3) * (3 - 2*math.log(4)) * (p1**-1) )
    b2 = mpf( (N2**2) * ((N2-1)**3) * (3 - 2*math.log(4)) * (p2**-1) )
    # define coefficient bounds of the quadratic taylor approximation
    # note that c1a <= c2a
    a1 = mpf( 0.5 )
    a2 = mpf( 0.5 * (N2 - 1) )
    # get the solution
    f = SumDefinitionMP(a1, b1, r1, a2, b2, r2)
    try:
        y0 = mp.findroot(
                f,
                [0, 1],
                solver=mp.calculus.optimization.Anderson,
                maxsteps=1000)
    except ValueError, e:
        return None
    print 'y0 =', y0
    #y0, r = scipy.optimize.brentq(f, 0.0, 1.0, full_output=True)
    #print f(0)
    #print f(1)
    #print y0
    #print 'root:', r.root
    #print 'iterations:', r.iterations
    #print 'function_calls:', r.function_calls
    #print 'converged:', r.converged
    # y0 = exp(-2 r2 t)
    try:
        t0 = -mp.log(y0) / (2 * r2)
    except ValueError, e:
        return None
    return t0

class SumDefinition:
    """
    This is for the numerical solution of an equation with a single root.
    The equation arises when we want to find the intersection of
    two functions which are each almost exponential functions.
    In particular, each one is the sum of a slowly decaying exponential
    and a faster exponential of particular forms,
    and the functions have particular relations to each other.
    This function has been transformed here so that x = exp(-t).
    Find the value of x that satisfies this inequality:
    a1 x^(2 r1) - b1 x^(3 r1) = a2 x^(2 r2) + b2 x*(3 r2)
    where
    0 < r1 < r2
    0 < a1 <= a2
    0 < b1
    0 < b2
    0 < x < 1
    make the substitution
    y = x^(2 r1)
    a1 y - b1 y^(3/2) = a2 y^(r2 / r1) + b2 y^((3/2)(r2/r1))
    divide by positive y
    a1 - b1 y^(1/2) = a2 y^(r2/r1 - 1) + b2 y^((3/2)(r2/r1) - 1)
    this function evaluates this difference
    """
    def __init__(self, a1, b1, r1, a2, b2, r2):
        """
        Initialize the constants.
        """
        error_message = 'bounds error'
        if not (0 < r1 < r2):
            raise ValueError(error_message)
        if not (0 < a1 <= a2):
            raise ValueError(error_message)
        if not (0 < b1):
            raise ValueError(error_message)
        if not (0 < b2):
            raise ValueError(error_message)
        self.a1 = a1
        self.b1 = b1
        self.r1 = r1
        self.a2 = a2
        self.b2 = b2
        self.r2 = r2
    def __call__(self, y):
        """
        This function should be used by the root finder.
        If the brentq root finder is used then the derivative does not help.
        """
        r = self.r2 / self.r1
        lhs = self.a1 - self.b1 * (y**0.5)
        rhs = self.a2 * (y**(r - 1)) + self.b2 * (y**(1.5*r - 1))
        return rhs - lhs

class SumDefinitionMP:
    """
    Multiple precision.
    """
    def __init__(self, a1, b1, r1, a2, b2, r2):
        """
        Initialize the constants.
        """
        error_message = 'bounds error'
        if not (0 < r1 < r2):
            raise ValueError(error_message)
        if not (0 < a1 <= a2):
            raise ValueError(error_message)
        if not (0 < b1):
            raise ValueError(error_message)
        if not (0 < b2):
            raise ValueError(error_message)
        self.a1 = a1
        self.b1 = b1
        self.r1 = r1
        self.a2 = a2
        self.b2 = b2
        self.r2 = r2
    def __call__(self, y):
        """
        This function should be used by the root finder.
        If the brentq root finder is used then the derivative does not help.
        """
        r = self.r2 / self.r1
        lhs = self.a1 - self.b1 * mp.root(y, 2)
        rhs = self.a2 * mp.power(y, r - 1) + self.b2 * mp.power(y, 1.5*r - 1)
        return rhs - lhs

def get_M(p):
    # This is the coefficient of a uniform bound.
    # M = 2 * mi_taylor_h(-p/2, p) / p
    # Now algebraically find this value in terms of p.
    # this is apparently 2 * (3/2 - log(4)) / p
    # which is (3 - 2 log(4)) / p
    return (3 - 2 * math.log(4)) / p

def mi_taylor_h(x, y):
    """
    Compute the Taylor expansion h function for mutual information.
    @param x: becomes small over time
    @param y: constant over time and depends on the stationary distribution
    """
    # compute the exact contribution to mutual information
    mi_full = y * (y + x) * math.log(1 + x / y)
    # compute the 2nd order Taylor approximate contribution
    mi_poly = x*y + 0.5 * x * x
    # compute the h function of the 2nd order approximation
    h = (mi_full - mi_poly) / (x * x)
    return h

def mi_taylor_h_grad(x, y):
    w = (2*y + x) * math.log(1 + x/y) - 2*x
    d_dx = - w * (x**-3) * y
    d_dy = w * (x**-2)
    return np.array([d_dx, d_dy])

def mi_cubic_abs(coeff, lam, t):
    """
    This is a pure exponential function.
    It is an upper bound on
    the absolute value of the sum of the O(x^3) Taylor error terms.
    The coefficient should be
    something like N^2 (N-1)^3 (3 - 2 log 4) (1/p).
    This is probably not going to work unless we switch to logarithms.
    But the value could be negative.
    @param coeff: the coefficient of the exponential function
    @param lam: the eigenvalue
    @param t: time
    """
    return coeff * math.exp(3 * lam * t)


#############################################
# Define functions for numerical optimization.
# The _mv suffix means multivariate.

def mi_taylor_h_mv(X):
    return mi_taylor_h(*X)

def mi_taylor_h_neg_mv(X):
    return -mi_taylor_h(*X)

def mi_taylor_h_grad_mv(X):
    return mi_taylor_h_grad(*X)

def mi_taylor_h_neg_grad_mv(X):
    return -mi_taylor_h_grad(*X)


def get_H_bound(x_abs_bound, y_low, y_high):
    """
    This explicitly computes the absolute bound of a particular function.
    The function is bivariate and the bounds of interest
    are over a given rectangular domain.
    The rectangular domain is subject to some constraints.
    @param x_abs_bound: bound on the absolute value of x
    @param y_low: lower bound on y
    @param y_high: upper bound on y
    @return: an abs bound on the function over the domain
    """
    if not (0 < y_low < y_high < 1):
        raise ValueError('invalid bounds on y')
    if not (0 < x_abs_bound < y_low):
        raise ValueError('invalid abs bound on x')
    x_low = -x_abs_bound
    x_high = x_abs_bound
    # Use a numerical minimizer as a sanity check.
    # Use an arbitrary point as the initial guess.
    t = .72
    # Minimize h.
    guess = np.array([
            (1 - t)*x_low + t*x_high,
            (1 - t)*y_low + t*y_high])
    result = scipy.optimize.fmin_l_bfgs_b(
            mi_taylor_h_mv,
            guess,
            fprime = mi_taylor_h_grad_mv,
            bounds = np.array([
                [x_low, x_high],
                [y_low, y_high]]))
    (x_argmin, y_argmin), h_min, info = result
    # Maximize h.
    guess = np.array([
            (1 - t)*x_low + t*x_high,
            (1 - t)*y_low + t*y_high])
    result = scipy.optimize.fmin_l_bfgs_b(
            mi_taylor_h_neg_mv,
            guess,
            fprime = mi_taylor_h_neg_grad_mv,
            bounds = np.array([
                [x_low, x_high],
                [y_low, y_high]]))
    (x_argmax, y_argmax), h_min_neg, info = result
    h_max = -h_min_neg
    """
    # report some values
    print 'x abs bound:'
    print x_abs_bound
    print
    print 'y lower bound:'
    print y_low
    print
    print 'y upper bound:'
    print y_high
    print
    print 'x that minimizes h:'
    print x_argmin
    print
    print 'y that minimizes h:'
    print y_argmin
    print
    print 'min h:'
    print h_min
    print
    print 'x that maximizes h:'
    print x_argmax
    print
    print 'y that maximizes h:'
    print y_argmax
    print
    print 'max h:'
    print h_max
    print
    """
    # check some properties of the numerical optimization
    if not (h_max > -h_min):
        msg = 'expected the max to have greater magnitude than the min'
        raise ValueError(msg)
    if not np.allclose(x_argmax, x_low):
        msg = 'expected the x argmax to be at the minimum allowed x'
        raise ValueError(msg)
    if not np.allclose(y_argmax, y_low):
        msg = 'expected the y argmax to be at the minimum allowed y'
        raise ValueError(msg)
    # return the max absolute value of h on the domain
    return h_max

class TestMe(unittest.TestCase):

    def test_h_bound(self):
        y_low = 0.1
        y_high = 0.9
        x_abs_bound = 0.09
        h = get_H_bound(x_abs_bound, y_low, y_high)


if __name__ == '__main__':
    unittest.main()
