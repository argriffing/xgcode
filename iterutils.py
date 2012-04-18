"""Iteration utilities.
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

def powerset(iterable):
    """
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    This is an itertools recipe from docs.python.org.
    """
    s = list(iterable)
    n = len(s)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(n+1))

def unique_everseen(iterable, key=None):
    """
    List unique elements, preserving order. Remember all elements ever seen.
    This is an itertools recipe from docs.python.org.
    """
    # unique_everseen('AAAABBBCCDAABBB') --> A B C D
    # unique_everseen('ABBCcAD', str.lower) --> A B C D
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in itertools.ifilterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def chopped(sequence, size):
    """
    Yields regular sized chunks of a sequence, but the last one may be ragged.
    Notice that strings are treated differently than other iterables.
    TODO replace this with an itertools recipe
    """
    if size <= 0:
        raise ValueError('not enough size')
    if len(sequence) == 1:
        yield sequence[0]
        return
    chunk = []
    for item in sequence:
        chunk.append(item)
        if len(chunk) == size:
            if type(sequence) in (str, unicode):
                yield ''.join(chunk)
            else:
                yield tuple(chunk)
            chunk = []
    if chunk:
        if type(sequence) in (str, unicode):
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

def ragged_grouper(seq, size):
    """
    This is a custom ragged variant of the grouper recipe.
    Unlike the chopper recipe this works with iterables.
    @param seq: a sequence with a length
    @param size: the size of the chunks
    """
    short_arr = []
    long_arr = []
    for v in seq:
        short_arr.append(v)
        if len(short_arr) == size:
            long_arr.append(tuple(short_arr))
            short_arr = []
    if short_arr:
        long_arr.append(tuple(short_arr))
    return long_arr

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

    def test_ragged_grouper(self):
        mygen = (a for a in range(10))
        observed = tuple(ragged_grouper(mygen, 3))
        expected = ((0,1,2),(3,4,5),(6,7,8),(9,))
        self.assertEqual(observed, expected)

    def test_rle(self):
        seq = (1, 1, 1, 2, 2, 3)
        result = tuple(rle(seq))
        self.assertEquals(result, ((1, 3), (2, 2), (3, 1)))



if __name__ == '__main__':
    unittest.main()
