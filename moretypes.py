"""
Simple types for argparse.
"""

def whole_number(x):
    try:
        x = int(x)
    except ValueError:
        raise TypeError
    if x < 1:
        raise TypeError
    return x

def positive_float(x):
    try:
        x = float(x)
    except ValueError:
        raise TypeError
    if x <= 0:
        raise TypeError
    return x

def int_ge_2(x):
    x = whole_number(x)
    if x < 2:
        raise TypeError
    return x
