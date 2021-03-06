"""
Constrained simple types for argparse.
"""

def _int(x):
    try:
        x = int(x)
    except ValueError:
        raise TypeError('expected an integer')
    return x

def _float(x):
    try:
        x = float(x)
    except ValueError:
        raise TypeError('expected a float')
    return x

def positive_integer(x):
    x = _int(x)
    if x < 1:
        raise TypeError('expected a positive integer')
    return x

def positive_float(x):
    x = _float(x)
    if x <= 0:
        raise TypeError('expected a positive float')
    return x

def nonneg_int(x):
    x = _int(x)
    if x < 0:
        raise TypeError('expected a non-negative integer')

def pos_int(x):
    x = _int(x)
    if x < 1:
        raise TypeError('expected a positive integer')
    return x

def nonneg_float(x):
    x = _float(x)
    if x < 0:
        raise TypeError('expected a non-negative float')

def pos_float(x):
    x = _float(x)
    if x <= 0:
        raise TypeError('expected a positive float')
    return x

class int_ge:
    def __init__(self, k):
        self.k = k
    def __call__(self, x):
        x = _int(x)
        if not (x >= self.k):
            raise TypeError(
                    'expected an integer '
                    'greater than or equal to %d' % self.k)
        return x

def whitespace_separated_sequence(s):
    seq = s.split()
    return tuple(seq)

