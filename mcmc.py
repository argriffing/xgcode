"""
Summarize MCMC samples of a statistic.

This uses formulas from BEAST Tracer by Andrew Rambaut et al.
Also it uses transliterations by Joseph Heled of the Tracer code.
"""

import math
import unittest

import numpy as np


g_max_lag = 2000


def get_hpd_interval(proportion, values):
    """
    Get the shortest interval that contains the right proportion.
    @param proportion: proportion of elements in the interval
    @param values: the sequence of posterior samples
    @return: low, high
    """
    d = sorted(values)
    ndata = len(d)
    nin = int(round(proportion * ndata))
    if nin < 2:
        raise ValueError('not enough samples')
    # k is the index of the lower endpoint
    # r is the width of the interval required to cover nin elements
    best_k = 0
    best_r = d[nin - 1] - d[0]
    for k in range(len(d) - (nin - 1)):
        r = d[k + nin - 1] - d[k]
        if r < best_r:
            best_r = r
            best_k = k
    return d[best_k], d[best_k + nin - 1]


class Correlation:
    def __init__(self):
        # mean
        self.mean = None
        # standard error of mean
        self.stdErrorOfMean = None
        # auto correlation time
        self.ACT = None
        # effective sample size
        self.ESS = None
        # standard deviation of autocorrelation time
        self.stdErrOfACT = None
    def analyze(self, values, stepsize=1):
        data = np.array(values)
        self.mean = np.mean(data)
        nsamples = len(data)
        nlags = min(nsamples-1, g_max_lag)
        gammas = [0] * nlags
        varstat = 0
        normalized = data - self.mean
        for lag in range(nlags):
            v1 = normalized[:nsamples - lag]
            v2 = normalized[lag:]
            v = v1 * v2
            gammas[lag] = np.mean(v)
            if lag == 0:
                varstat = gammas[0]
            elif lag % 2 == 0:
                # fancy stopping criterion
                if gammas[lag - 1] + gammas[lag] > 0:
                    varstat += 2 * (gammas[lag - 1] + gammas[lag])
                else:
                    nlags = lag
                    break
        self.stdErrorOfMean = math.sqrt(varstat / nsamples)
        self.ACT = stepsize * varstat / gammas[0]
        self.ESS = stepsize * nsamples / self.ACT
        coeff = 2 * (varstat / gammas[0]) * stepsize
        self.stdErrOfACT = coeff * math.sqrt(4 * (nlags + 1) / float(nsamples))


class TestStats(unittest.TestCase):

    def test_hpd(self):
        values = (
                0.05, 0.4, 0.25, 2.44, 19.8, 7.26, 1, 0.3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3,
                0.1, 0.3, 0.2, 3.4, 9.8, 7.6, 1.1, 0.3)
        interval = get_hpd_interval(0.95, values)
        self.assertEqual(interval, (0.1, 3))


if __name__ == '__main__':
    unittest.main()

