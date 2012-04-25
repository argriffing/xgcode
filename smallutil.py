"""
A small set of functions.
"""

import itertools

def pairwise(iterable):
    """
    Yield pairs of neighbors.
    This is directly from the python itertools documentation.
    For example,
    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itertools.tee(iterable)
    next(b, None)
    return itertools.izip(a, b)

def stripped_lines(lines):
    """
    This function yields nonempty stripped lines.
    """
    for line in lines:
        line = line.strip()
        if line:
            yield line

def get_stripped_lines(lines):
    """
    @return: a list of nonempty stripped lines
    """
    return list(stripped_lines(lines))

