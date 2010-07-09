"""
Run an algorithm until one of three conditions happens.
First, the user could manually break out of the loop with ctrl-c.
Second, the loop could hit a time limit.
Third, the loop could hit a limit on the number of iterations.
"""

import time


class ComboBreaker(Exception): pass


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


def combo_breaker(state, fn, nseconds=None, niterations=None):
    """
    Raise a ComboBreaker exception containing the final state.
    @param state: a state
    @param fn: functionally update state using this function
    @param nseconds: a time limit
    @param niterations: a limit on the number of iterations
    """
    tm = time.time()
    count = 0
    try:
        try:
            while True:
                if nseconds is not None:
                    if time.time() - tm >= nseconds:
                        raise ComboBreaker('time limit')
                if niterations is not None:
                    if count >= niterations:
                        raise ComboBreaker('iteration limit')
                state = fn(state)
                count += 1
        except KeyboardInterrupt:
            raise ComboBreaker('ctrl-c')
    except ComboBreaker as e:
        return RunInfo(state, count, time.time() - tm, str(e))

