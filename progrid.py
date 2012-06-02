"""
This is a module for progressive gridding.
The idea is to yield points in an interval in a way that
can be stopped at any time.
For a 1d interval the endpoints will be yielded first.
An example application is plotting points.
"""

import random
import unittest
import heapq

import cbreaker


def gen_binary(low, high):
    yield low
    yield high
    k = 1
    while True:
        n = 2**k
        incr = (high - low) / float(n)
        for i in range(1, n):
            if i % 2:
                yield low + i * incr
        k += 1

def gen_poisson(low, high):
    yield low
    yield high
    while True:
        yield random.uniform(low, high)

def gen_overdispersed(low, high):
    """
    This is a generator that samples overdispersed events in an interval.
    It requires a lot of memory.
    """
    # Create the min heap.
    # The triples are the neg length, the low, and the high values.
    # Popping from the queue will return the longest interval.
    # The queue grows linearly with the number of sampled events.
    yield low
    yield high
    q = [(low-high, low, high)]
    while True:
        dummy, a, b = heapq.heappop(q)
        mid = random.uniform(a, b)
        heapq.heappush(q, (a-mid, a, mid))
        heapq.heappush(q, (mid-b, mid, b))
        yield mid


class TestProGrid(unittest.TestCase):

    def test_binary(self):
        observed = list(cbreaker.throttled(gen_binary(0, 1), ncount=9))
        expected = [0.0, 1.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875]
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()

