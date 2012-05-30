"""
This module has a stable implementation of the Bernoulli generating function.

This generating function happens to be related to the construction
of a mutation-selection balance rate matrix
given allele selection coefficients and a mutation rate matrix.
The exponential generating function of the Berouli numbers will be denoted bgf.
http://mathworld.wolfram.com/BernoulliNumber.html
"""

import unittest
import math

def bgf(x):
    """
    This function computes x / (e^x - 1) in a numerically stable way.
    @param x: an unrestricted floating point number
    @return: a positive float
    """
    if abs(x) < 0.01:
        return bgf_near_zero(x)
    else:
        return bgf_naive_expm1(x)

def bgf_naive(x):
    return x / (math.exp(x) - 1)

def bgf_naive_expm1(x):
    """
    The expm1 function is a C standard so it should be thinly wrapped by Python.
    Its purpose is accuracy near zero but that feature is not used here,
    because the ratio is badly behaved near zero anyway.
    """
    try:
        denominator = math.expm1(x)
    except OverflowError as e:
        return 0
    else:
        return x / denominator

def bgf_near_zero(x):
    """
    This is the series expansion at x=0.
    The error term is O(x^8).
    It is accurate beyond all reasonable expectation when abs(x) < 0.01.
    This series expansion and its error characteristics
    are from a computer algebra system on the internet.
    The Horner-like evaluation is my own factorization.
    This could probably benefit from some kind of C extension
    or Cython treatment.
    """
    return 1 - (x/2)*(1 - (x/6)*(1 - (x*x/60)*(1 - (x*x/42))))


class TestBgf(unittest.TestCase):

    def test_bgf_stable(self):
        eps = 1e-8
        for sign in (-1, 1):
            for x in (0.0001, 0.1, 2.0, 333.0):
                y_bgf = bgf(x)
                y_bgf_naive = bgf_naive(x)
                self.assertTrue(
                        abs(y_bgf_naive - y_bgf) < eps,
                        (y_bgf, y_bgf_naive))

    def test_bgf_unstable(self):
        eps = 1e-8
        x = 1e-10
        y_bgf = bgf(x)
        y_bgf_naive = bgf_naive(x)
        self.assertFalse(
                abs(y_bgf_naive - y_bgf) < eps,
                (y_bgf, y_bgf_naive))

    def test_bgf_huge(self):
        self.assertEqual(bgf(800), 0)


if __name__ == '__main__':
    unittest.main()

