"""
Run an algorithm until one of four conditions happens.
First, the user could manually break out of the loop with ctrl-c.
Second, the loop could hit a time limit.
Third, the loop could hit a limit on the number of iterations.
Fourth, the algorithm could terminate successfully.
This framework is intended to support open ended computations.
For example, the computation could look for a counterexample to a conjecture.
Or it could compute a statistic on some arbitrarily large set of sampled data.
"""

import time
import unittest


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
        chunks = [
            'termination condition: %s' % self.stop_msg,
            'completed iterations: %s' % self.iterations,
            'elapsed seconds: %s' % self.seconds,
            self.state.get_response()]
        return '\n'.join(chunks)


def combo_breaker(states, nseconds=None, niterations=None):
    """
    Return a RunInfo object containing the final state.
    @param states: a state iterator, probably a generator object
    @param nseconds: a time limit
    @param niterations: a limit on the number of iterations
    """
    tm = time.time()
    count_m1 = -1
    state = None
    stop_msg = 'finished'
    try:
        for count_m1, state in enumerate(states):
            if nseconds is not None:
                if time.time() - tm >= nseconds:
                    stop_msg = 'time limit'
                    break
            if niterations is not None:
                if count_m1 + 1 >= niterations:
                    stop_msg = 'iteration limit'
                    break
    except KeyboardInterrupt:
        stop_msg = 'ctrl-c'
    return RunInfo(state, count_m1 + 1, time.time() - tm, stop_msg)


class TestComboBreaker(unittest.TestCase):

    def test_combo_breaker_finish(self):
        info = combo_breaker(xrange(3))
        self.assertEqual(info.iterations, 3)
        self.assertEqual(info.stop_msg, 'finished')

    def test_combo_breaker_niterations(self):
        info = combo_breaker(xrange(3), niterations=2)
        self.assertEqual(info.iterations, 2)
        self.assertEqual(info.stop_msg, 'iteration limit')


if __name__ == '__main__':
    unittest.main()
