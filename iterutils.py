"""
Iteration utilities.
The name of this module is semi-standard.
See the comment by Alex Martelli.
http://stackoverflow.com/questions/1639772
converting-a-single-ordered-list-in-python-to-a-dictionary-pythonically
Each function should use constant memory.
Some are standard itertools recipes.
Some should be converted to standard itertools recipes.
Some are miscellaneous custom functions.
"""

import itertools
import unittest


class Filler:
    """
    Generate values over a range of sequential positions.
    At some of the positions a value is available,
    but at other positions no value is available.
    """
    def __init__(self, low, high, default_value, truncate=False):
        """
        @param low: the position for which the first value is yielded
        @param high: the position for which the last value is yielded
        @param default_value: the value used for filler
        @param truncate: disregard positions which are too low or too high
        """
        self.low = low
        self.high = high
        self.default_value = default_value
        self.truncate = truncate
        self.prev = None

    def fill(self, position, value):
        """
        Yield an informative value and maybe some uninformative ones.
        This function should be called repeatedly,
        and with strictly increasing positions.
        @param position: an available position
        @param value: the value at the position
        """
        # If the position is outside the range, 
        # deal with it according to the truncation option.
        if not (self.low <= position <= self.high):
            if self.truncate:
                return
            else:
                msg_a = 'position %d ' % position
                msg_b = 'is outside [%d, %d]' % (self.low, self.high)
                raise ValueError(msg_a + msg_b)
        # check monotonicity
        if self.prev is not None:
            if position <= self.prev:
                raise ValueError('positions should monotonically increase')
        # fill between the previous position and the current position
        for i in xrange(self.get_ngap(position)):
            yield default_value
        # yield the value at the current position
        yield value
        self.prev = position

    def finish(self):
        nremaining = self.get_ngap(self.high) + 1
        for i in xrange(self.get_nremaining()):
            yield default_value
        self.prev = self.high

    def get_ngap(self, position):
        if self.prev is None:
            return position - self.low
        else:
            return (position - self.prev) - 1


def get_only(collection):
    """
    @param collection: a collection with only one item
    """
    if len(collection) != 1:
        raise ValueError('expected the collection to have a single element')
    for value in collection:
        return value

def read_backwards(fin, blocksize=4096):
    """
    Read a file line by line, backwards.
    http://code.activestate.com/recipes/439045/
    contrib: Peter Astrand, Raymond Hettinger
    @param fin: a file open for reading
    @param blocksize: read this many bytes at a time
    """
    buf = ""
    fin.seek(-1, 2)
    lastchar = fin.read(1)
    trailing_newline = (lastchar == "\n")
    while 1:
        newline_pos = buf.rfind("\n")
        pos = fin.tell()
        if newline_pos != -1:
            # found a newline
            line = buf[newline_pos+1:]
            buf = buf[:newline_pos]
            if pos or newline_pos or trailing_newline:
                line += "\n"
            yield line
        elif pos:
            # need to fill buffer
            toread = min(blocksize, pos)
            fin.seek(-toread, 1)
            buf = fin.read(toread) + buf
            fin.seek(-toread, 1)
            if pos == toread:
                buf = "\n" + buf
        else:
            # start of file
            return

def stripped_lines(lines):
    """
    This function yields nonempty stripped lines.
    """
    for line in lines:
        line = line.strip()
        if line:
            yield line

def dot_product(va, vb):
    return sum(a * b for a, b in itertools.izip(va, vb))

def product(numbers):
    x = 1
    for number in numbers:
        x *= number
    return x

def rle(sequence):
    """
    Yield (value, length) pairs.
    TODO replace this with an itertools recipe
    """
    first = True
    for v in sequence:
        if first:
            first = False
            value = v
            count = 1
        elif value == v:
            count += 1
        else:
            yield (value, count)
            value = v
            count = 1
    yield (value, count)

def chopped(sequence, size):
    """
    Yields regular sized chunks of a sequence, but the last one may be ragged.
    Notice that strings are treated differently than other iterables.
    TODO replace this with an itertools recipe
    """
    assert size > 0
    if len(sequence) == 1:
        yield sequence[0]
        return
    chunk = []
    for item in sequence:
        chunk.append(item)
        if len(chunk) == size:
            if type(sequence) is str:
                yield ''.join(chunk)
            else:
                yield tuple(chunk)
            chunk = []
    if chunk:
        if type(sequence) is str:
            yield ''.join(chunk)
        else:
            yield tuple(chunk)

def grouper(sequence, size, fillvalue=None):
    """
    This is like the chopped function but works for generators.
    This is directly from the python itertools documentation.
    For example,
    grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    @param sequence: an iterable
    @param size: the size of the chunks
    """
    args = [iter(sequence)] * size
    return itertools.izip_longest(fillvalue=fillvalue, *args)

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

def chopped_nonbreaking(sequence, size):
    """
    Like L{chopped} but forces consecutive identical elements to be in the same chunk.
    Chunks containing consecutive identical elements may be larger than the desired chunk size.
    The last chunk may be smaller than the desired chunk size.
    @param sequence: something iterable
    @param size: the desired size of a chunk
    """
    assert size > 0
    if len(sequence) == 1:
        yield sequence[0]
        return
    chunk = []
    for item in sequence:
        # if the chunk has at least the desired size
        # and the item is different from the last chunk item,
        # then start a new chunk.
        if len(chunk) >= size and chunk[-1] != item:
            if type(sequence) is str:
                yield ''.join(chunk)
            else:
                yield tuple(chunk)
            chunk = []
        chunk.append(item)
    if chunk:
        if type(sequence) is str:
            yield ''.join(chunk)
        else:
            yield tuple(chunk)


class TestIterutils(unittest.TestCase):

    def test_chopped_nonbreaking(self):
        seq = (1, 1, 1, 2, 2, 3)
        groups = tuple(chopped_nonbreaking(seq, 3))
        self.assertEquals(groups, ((1, 1, 1), (2, 2, 3)))
        groups = tuple(chopped_nonbreaking(seq, 2))
        self.assertEquals(groups, ((1, 1, 1), (2, 2), (3,)))

    def test_chopped(self):
        seq = (1, 1, 1, 2, 2, 3)
        groups = tuple(chopped(seq, 3))
        self.assertEquals(groups, ((1, 1, 1), (2, 2, 3)))
        groups = tuple(chopped(seq, 2))
        self.assertEquals(groups, ((1, 1), (1, 2), (2, 3)))

    def test_grouper(self):
        mygen = (a for a in range(10))
        observed = tuple(grouper(mygen, 3))
        expected = ((0,1,2),(3,4,5),(6,7,8),(9,None,None))
        self.assertEqual(observed, expected)

    def test_rle(self):
        seq = (1, 1, 1, 2, 2, 3)
        result = tuple(rle(seq))
        self.assertEquals(result, ((1, 3), (2, 2), (3, 1)))

    def test_dot_product(self):
        self.assertEquals(dot_product((1, 2, 3), (4, 5, 6)), 1*4 + 2*5 + 3*6)


if __name__ == '__main__':
    unittest.main()
