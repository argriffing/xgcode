"""
Validation functions to use with argparse.
"""

def positive_integer(s):
    i = int(s)
    if i < 1:
        raise TypeError('expected a positive integer')
    return i

def positive_float(s):
    i = float(s)
    if i <= 0:
        raise TypeError('expected a positive float')
    return i
