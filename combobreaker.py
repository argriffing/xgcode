"""
Run an algorithm until one of four conditions happens.
First, the user could manually break out of the loop with Ctrl-C.
Second, the loop could hit a time limit.
Third, the loop could hit a limit on the number of iterations.
Fourth, the algorithm could terminate successfully.
This framework is intended to support open ended computations.
For example, the computation could look for a counterexample to a conjecture.
Or it could compute a statistic on some arbitrarily large set of sampled data.
"""

import time
import unittest

STOP_FINISHED = 'finished'
STOP_ITERATION_LIMIT = 'iteration limit'
STOP_TIME_LIMIT = 'time limit'
STOP_INTERRUPT = 'Ctrl-C'


class RunInfo(object):
    """
    Provide information about the run.
    """
    def __init__(self, state, iterations, seconds, stop_msg):
        """
        @param state: user defined state
        @param iterations: the number of iterations completed
        @param seconds: the number of seconds elapsed
        @param stop_msg: the termination message
        """
        self.state = state
        self.iterations = iterations
        self.seconds = seconds
        self.stop_msg = stop_msg

    def get_response(self):
        if self.state is None:
            state_response = 'no state available'
        else:
            state_response = self.state.get_response()
        return '\n'.join([
            'termination condition: %s' % self.stop_msg,
            'completed iterations: %s' % self.iterations,
            'elapsed seconds: %s' % self.seconds,
            state_response])


def run(states, nseconds=None, niterations=None):
    """
    Return a RunInfo object containing the final state.
    @param states: a state iterator, probably a generator object
    @param nseconds: a time limit
    @param niterations: a limit on the number of iterations
    """
    tm = time.time()
    count_m1 = -1
    state = None
    stop_msg = STOP_FINISHED
    try:
        for count_m1, state in enumerate(states):
            if nseconds is not None:
                if time.time() - tm >= nseconds:
                    stop_msg = STOP_TIME_LIMIT
                    break
            if niterations is not None:
                if count_m1 + 1 >= niterations:
                    stop_msg = STOP_ITERATION_LIMIT
                    break
    except KeyboardInterrupt:
        stop_msg = STOP_INTERRUPT
    return RunInfo(state, count_m1 + 1, time.time() - tm, stop_msg)

class Collatz:
    """
    For testing.
    """
    def __init__(self, k):
        self.k = k
        self.seq = [k]
        self.high_water = [k]
    def __iter__(self):
        return self
    def next(self):
        if self.k == 1:
            raise StopIteration
        elif self.k % 2:
            self.k = 3*self.k + 1
        else:
            self.k /= 2
        self.seq.append(self.k)
        if self.k > self.high_water[-1]:
            self.high_water.append(self.k)
        return self
    def get_response(self):
        return '\n'.join([
            'sequence: ' + str(self.seq),
            'high water marks: ' + str(self.high_water)])


class TestComboBreaker(unittest.TestCase):

    def test_finish(self):
        info = run(xrange(3))
        self.assertEqual(info.iterations, 3)
        self.assertEqual(info.stop_msg, STOP_FINISHED)

    def test_niterations(self):
        info = run(xrange(3), niterations=2)
        self.assertEqual(info.iterations, 2)
        self.assertEqual(info.stop_msg, STOP_ITERATION_LIMIT)

    def test_collatz(self):
        info = run(Collatz(42), niterations=5)
        self.assertEqual(info.iterations, 5)
        self.assertEqual(info.stop_msg, STOP_ITERATION_LIMIT)

    def test_collatz_premature(self):
        info = run(Collatz(42))
        self.assertEqual(info.iterations, 8)
        self.assertEqual(info.stop_msg, STOP_FINISHED)
        self.assertEqual(len(info.get_response().splitlines()), 5)

    def test_collatz_degenerate(self):
        info = run(Collatz(1))
        self.assertEqual(info.iterations, 0)
        self.assertEqual(info.stop_msg, STOP_FINISHED)


if __name__ == '__main__':
    unittest.main()

