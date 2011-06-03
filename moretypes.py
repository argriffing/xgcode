"""
Simple types for argparse.
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

def whole_number(x):
    x = _int(x)
    if x < 1:
        raise TypeError
    return x

def positive_float(x):
    x = _float(x)
    if x <= 0:
        raise TypeError
    return x

class int_ge:
    def __init__(k):
        self.k = k
    def __call__(self, x):
        x = _int(x)
        if not (x >= self.k):
            msg = 'expected an integer greater than or equal to %d' % self.k
            raise TypeError(msg)
        return x
