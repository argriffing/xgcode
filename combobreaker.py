"""
Run an algorithm until one of four conditions happens.
First, the user could manually break out of the loop with Ctrl-C.
Second, the loop could hit a time limit.
Third, the loop could hit a limit on the number of iterations.
Fourth, the algorithm could terminate successfully.
This framework is intended to support open ended computations.
For example, the computation could look for a counterexample to a conjecture.
Or it could compute a statistic on some arbitrarily large set of sampled data.
Currently three primary interfaces are available.
These are run_callable, run_checker, and run.
The last of these is obsolescent.
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

    def __str__(self):
        lines = [
            'loop termination explanation: %s' % self.stop_msg,
            'completed iterations: %s' % self.iterations,
            'elapsed seconds: %s' % self.seconds,
            str(self.state)]
        return '\n'.join(lines)


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

def run_checker(
        checker, values,
        nseconds=None, niterations=None):
    """
    Return a RunInfo object containing the final state.
    @param states: a state iterator, probably a generator object
    @param nseconds: a time limit
    @param niterations: a limit on the number of iterations
    """
    tm = time.time()
    stop_msg = STOP_FINISHED
    nchecked = 0
    try:
        for value in values:
            found = checker(value)
            nchecked += 1
            if found:
                break
            if nseconds is not None:
                if time.time() - tm >= nseconds:
                    stop_msg = STOP_TIME_LIMIT
                    break
            if niterations is not None:
                if nchecked >= niterations:
                    stop_msg = STOP_ITERATION_LIMIT
                    break
    except KeyboardInterrupt:
        stop_msg = STOP_INTERRUPT
    return RunInfo(checker, nchecked, time.time() - tm, stop_msg)

def run_callable(f, nseconds=None, niterations=None):
    """
    Return a RunInfo object containing the final state.
    The callable returns true if the search has finished.
    @param states: a state iterator, probably a generator object
    @param nseconds: a time limit
    @param niterations: a limit on the number of iterations
    """
    tm = time.time()
    stop_msg = STOP_FINISHED
    nchecked = 0
    try:
        while True:
            found = f()
            nchecked += 1
            if found:
                break
            if nseconds is not None:
                if time.time() - tm >= nseconds:
                    stop_msg = STOP_TIME_LIMIT
                    break
            if niterations is not None:
                if nchecked >= niterations:
                    stop_msg = STOP_ITERATION_LIMIT
                    break
    except KeyboardInterrupt:
        stop_msg = STOP_INTERRUPT
    return RunInfo(f, nchecked, time.time() - tm, stop_msg)

def collatz(k):
    """
    This is for testing, and it runs forever.
    @param k: the starting value which should be a positive integer
    """
    while True:
        yield k
        if k % 2:
            k = 3*k + 1
        else:
            k /= 2

def collatz_checker(k):
    """
    This is for testing.
    @param k: a value to check
    @return: True when k satisfies the acceptance conditions
    """
    return k == 1

class CollatzTrackingChecker:
    """
    For testing, this looks for collatz hitting 1 for a given start value.
    This is essentially a stateful version of collatz_checker.
    """
    def __init__(self):
        self.seq = []
        self.high_water = []
    def __call__(self, k):
        """
        @param k: a value to check
        @return: True when k satisfies the acceptance conditions
        """
        self.seq.append(k)
        if not self.high_water or k > self.high_water[-1]:
            self.high_water.append(k)
        return k == 1
    def __str__(self):
        return '\n'.join([
            'sequence: ' + str(self.seq),
            'high water marks: ' + str(self.high_water)])

class CollatzCallable:
    """
    For testing.
    """
    def __init__(self, k):
        self.initial_k = k
        self.k = None
        self.seq = []
        self.high_water = []
    def __call__(self):
        if self.k is None:
            self.k = self.initial_k
        elif self.k % 2:
            self.k = 3*self.k + 1
        else:
            self.k /= 2
        self.seq.append(self.k)
        if not self.high_water or self.k > self.high_water[-1]:
            self.high_water.append(self.k)
        return self.k == 1
    def __str__(self):
        return '\n'.join([
            'sequence: ' + str(self.seq),
            'high water marks: ' + str(self.high_water)])

class Collatz:
    """
    For testing.
    This is for an obsolescent interface.
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
    def __str__(self):
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
        info = run(Collatz(42))
        self.assertEqual(info.iterations, 8)
        self.assertEqual(info.stop_msg, STOP_FINISHED)
        self.assertEqual(len(str(info).splitlines()), 5)

    def test_collatz_premature(self):
        info = run(Collatz(42), niterations=5)
        self.assertEqual(info.iterations, 5)
        self.assertEqual(info.stop_msg, STOP_ITERATION_LIMIT)

    def test_collatz_degenerate(self):
        info = run(Collatz(1))
        self.assertEqual(info.iterations, 0)
        self.assertEqual(info.stop_msg, STOP_FINISHED)
        self.assertEqual(len(str(info).splitlines()), 4)

    def test_collatz_nontracking_checker(self):
        info = run_checker(collatz_checker, collatz(42))
        response = str(info)
        self.assertEqual(info.iterations, 9)
        self.assertEqual(info.stop_msg, STOP_FINISHED)
        self.assertEqual(len(response.splitlines()), 4)

    def test_collatz_tracking_checker(self):
        info = run_checker(CollatzTrackingChecker(), collatz(42))
        self.assertEqual(info.iterations, 9)
        self.assertEqual(info.stop_msg, STOP_FINISHED)
        self.assertEqual(len(str(info).splitlines()), 5)

    def test_collatz_callable(self):
        info = run_callable(CollatzCallable(42))
        self.assertEqual(info.iterations, 9)
        self.assertEqual(info.stop_msg, STOP_FINISHED)
        self.assertEqual(len(str(info).splitlines()), 5)


if __name__ == '__main__':
    unittest.main()

