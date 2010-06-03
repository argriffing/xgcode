"""
Add stuff specific to this PCA fungus project.
"""

import numpy as np

import Util

class WordError(Exception): pass

class Word(object):
    def __init__(self, line):
        # store the original line
        self.line = line
        # extract the name and the binary vector
        elements = line.split()
        self.name = elements[0]
        for x in elements[1:]:
            if x not in ('0', '1'):
                msg = 'expected 0 or 1 but got ' + x
                raise WordError(msg)
        self.v = np.array([int(x) for x in elements[1:]])

def validate_words(words):
    """
    Make sure the binary vectors are the same lengths.
    @param words: word objects
    """
    lengths = set(len(w.v) for w in words)
    if len(lengths) != 1:
        msg = 'each binary vector should be the same length'
        raise WordError(msg)

def get_words(lines):
    """
    @param lines: lines of a .hud file
    @return: a list of validated words
    """
    lines = Util.get_stripped_lines(lines)
    words = [Word(line) for line in lines]
    validate_words(words)
    return words
