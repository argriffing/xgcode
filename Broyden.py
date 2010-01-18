import unittest
import numpy
import math
import sys

import scipy.optimize


def easy_function(X):
    return tuple(x*x - 1 for x in X)

class MyTest(unittest.TestCase):

    def helper(self, solve):
        guess = (0, 0, 0)
        iterations = 20
        result = solve(easy_function, guess, iterations)
        for x in result:
            self.assertAlmostEqual(x*x, 1.0)

    def test_current(self):
        "The current version of broyden2 fails."
        self.helper(scipy.optimize.nonlin.broyden2)

    def test_modified(self):
        "A modified version of broyden2 does not fail so hard."
        def norm(v):
            return math.sqrt(numpy.sum((numpy.array(v)**2).flat))
        def myF(F,xm):
            return numpy.matrix(F(tuple(xm.flat))).T
        def broyden2(F, xin, iter=10, alpha=0.4, verbose = False):
            xm=numpy.matrix(xin).T
            Fxm=myF(F,xm)
            Gm=-alpha*numpy.matrix(numpy.identity(len(xin)))
            for n in range(iter):
                deltaxm=-Gm*Fxm
                xm=xm+deltaxm
                Fxm1=myF(F,xm)
                # begin modification
                deltaFxm=Fxm1-Fxm
                denominator = norm(deltaFxm)**2
                if not denominator:
                    break
                Fxm=Fxm1
                Gm=Gm+(deltaxm-Gm*deltaFxm)*deltaFxm.T/denominator
                # end modification
                if verbose:
                    print "%d:  |F(x)|=%.3f"%(n+1, norm(Fxm))
            return xm.flat
        self.helper(broyden2)

if __name__ == '__main__':
    unittest.main()

