"""
Provide some simple functors that sample branch lengths.
Each functor is callable, returning a sampled branch length.
Each functor describes itself with str.
"""

import math
import random

class ShortAscii:
    """
    Sample branch lengths that do not take much space when written.
    """

    def __init__(self, allow_integers, allow_reciprocals):
        """
        @param allow_integers: allow single digit positive integers
        @param allow_reciprocals: allow reciprocals of single digit positive integers
        """
        self.pool = [1.0]
        if allow_integers:
            self.pool.extend([float(x) for x in range(2, 10)])
        if allow_reciprocals:
            self.pool.extend([1.0 / x for x in range(2, 10)])

    def __str__(self):
        return 'branch length are single digit positive integers or their reciprocals'

    def __call__(self):
        """
        @return: a float that is a positive single digit integer or its reciprocal
        """
        return random.choice(self.pool)

class Exponential:

    def __str__(self):
        return 'branch lengths are exponentially distributed with mean 0.1'

    def __call__(self):
        """
        @return: a branch length drawn from an exponential distribution
        """
        mean = 0.1
        return random.expovariate(1/mean)

class Pachter:
    """
    This branch length of 0.1 is used in the "Why neighbor joining works" paper.
    Lior Pacher is a co-author of this paper.
    """

    def __str__(self):
        return 'each branch has length 0.1'

    def __call__(self):
        """
        @return: a branch length drawn from a constant distribution
        """
        return 0.1

class Uniform:

    def __init__(self, low, high):
        """
        @param low: the low bound for uniform sampling
        @param high: the high bound for uniform sampling
        """
        self.low = low
        self.high = high

    def __str__(self):
        return 'branch lengths are uniformly distributed in [%s, %s]' % (self.low, self.high)

    def __call__(self):
        """
        @return: a branch length drawn from a uniform distribution
        """
        return random.uniform(self.low, self.high)

class UniformA(Uniform):
    def __init__(self):
        Uniform.__init__(self, 0, 1)

class UniformB(Uniform):
    def __init__(self):
        Uniform.__init__(self, 0.05, 0.15)
