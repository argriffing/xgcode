"""
An updated module to eventually replace the combobreaker module.
This is about running loops that terminate under any of various conditions,
for example on ctrl-c or timeout or iteration limit,
or when a sentinel value is reached or the search space is exhausted.
Instantiate a 'Throttled' object if you care why the termination stopped.
Otherwise use the 'throttled' wrapper.
If you need to accumulate a summary over the course of the iterations,
then use the more complicated 'combobreaker' module instead.
"""

from StringIO import StringIO
import time
import unittest

STOP_FINISHED = 'finished'
STOP_ITERATION_LIMIT = 'iteration limit'
STOP_TIME_LIMIT = 'time limit'
STOP_INTERRUPT = 'Ctrl-C'

class Stop(Exception):
    def __init__(self, condition):
        self.condition = condition

class Throttled:
    """
    Use this if you care why the iteration stopped.
    """
    def __init__(self):
        self.ncompleted = 0
        self.nseconds = None
        self.condition = None
    def run(self, states, nseconds=None, ncount=None):
        """
        @param states: an iterable
        @param nseconds: a time limit
        @param ncount: a limit on the number of yielded things
        """
        self.t0 = time.time()
        try:
            for state in throttled_raises(states, nseconds, ncount):
                yield state
                self.ncompleted += 1
        except Stop as e:
            self.condition = e.condition
            self.nseconds = time.time() - self.t0
    def __str__(self):
        out = StringIO()
        print >> out, 'loop termination explanation:', self.condition
        print >> out, 'completed iterations:', self.ncompleted
        print >> out, 'elapsed time:', self.nseconds
        return out.getvalue().rstrip()


def throttled(states, nseconds=None, ncount=None):
    """
    Use this if you do not care why the iteration stopped.
    @param states: an iterable
    @param nseconds: a time limit
    @param ncount: a limit on the number of yielded things
    """
    try:
        for x in throttled_raises(states, nseconds, ncount):
            yield x
    except Stop as e:
        pass

def throttled_raises(states, nseconds=None, ncount=None):
    """
    @param states: an iterable
    @param nseconds: a time limit
    @param ncount: a limit on the number of yielded things
    """
    it = iter(states)
    t0 = time.time()
    ncompleted = 0
    try:
        while True:
            if nseconds is not None:
                if time.time() - t0 >= nseconds:
                    raise Stop(STOP_TIME_LIMIT)
            if ncount is not None:
                if ncompleted >= ncount:
                    raise Stop(STOP_ITERATION_LIMIT)
            try:
                yield next(it)
            except StopIteration:
                raise Stop(STOP_FINISHED)
            ncompleted += 1
    except KeyboardInterrupt:
        raise Stop(STOP_INTERRUPT)


class TestCbreaker(unittest.TestCase):

    def test_throttled(self):
        self.assertEqual(
                list(throttled(range(20), ncount=10)),
                range(10))


if __name__ == '__main__':
    unittest.main()
