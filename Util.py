import random
import unittest

import iterutils

def get_first(elements):
    for element in elements:
        return element

def get_paragraphs(raw_lines):
    """
    Leading and trailing whitespace of each line is removed.
    Paragraphs are lists of lines that are separated by blank lines.
    @param raw_lines: raw lines of input
    @return: a list of stripped line lists
    """
    paragraphs = []
    paragraph = []
    for line in raw_lines:
        stripped_line = line.strip()
        if stripped_line:
            paragraph.append(stripped_line)
        else:
            if paragraph:
                paragraphs.append(paragraph)
            paragraph = []
    if paragraph:
        paragraphs.append(paragraph)
    return paragraphs

def get_stripped_lines(raw_lines):
    """
    @return: a list of nonempty stripped lines
    """
    return list(iterutils.stripped_lines(raw_lines))

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke.
    This is not by me.
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

def weights_to_distribution(weights):
    for weight in weights:
        assert weight >= 0
    weight_list = list(weights)
    total = sum(weight_list)
    assert total > 0
    return [weight / total for weight in weight_list]

def flattened_nonrecursive(lists):
    """
    @param lists: a list of lists
    @return: a list
    """
    arr = []
    for v in lists:
        arr.extend(v)
    return arr

def hamming_distance(first, second):
    return sum(1 for a, b in zip(first, second) if a != b)

def bresenham_line(x_count, y_count):
    """
    Yield a y value for each x.
    This function can be used to chop sequences.
    @param x_count: the number of x pixels
    @param y_count: the number of y pixels
    """
    dx = int(x_count - 1)
    dy = int(y_count - 1)
    assert dx == x_count - 1
    assert dy == y_count - 1
    assert dx > 0
    assert dy > 0
    assert dx >= dy
    error = 0.0
    deltaerr = dy / float(dx)
    y = 0
    for x in range(x_count):
        yield y
        error += deltaerr
        if abs(error) >= 0.5:
            y += 1
            error -= 1.0

def chopped_bresenham(sequence, nchunks):
    """
    Yield chunks of a sequence.
    The lengths of the yielded chunks differ from each other by at most one.
    @param sequence: can be iterated and its length can be found
    @param nchunks: the number of yielded chunks
    """
    assert len(sequence) >= nchunks
    chunk = []
    current_category = None
    for v, category in zip(sequence, bresenham_line(len(sequence), nchunks)):
        if category == current_category:
            chunk.append(v)
        else:
            if current_category is not None:
                if type(sequence) is str:
                    yield ''.join(chunk)
                else:
                    yield tuple(chunk)
            current_category = category
            chunk = [v]
    if chunk:
        if type(sequence) is str:
            yield ''.join(chunk)
        else:
            yield tuple(chunk)

def random_weighted_int(weights):
    """
    @param weights: an ordered sequence of nonnegative weights summing to one
    """
    x = random.random()
    accum = 0
    for i, w in enumerate(weights):
        accum += w
        if x < accum:
            return i

def weighted_choice(weight_state_pairs):
    if not weight_state_pairs:
        raise ValueError('no choices available')
    if len(weight_state_pairs) == 1:
        weight, state = weight_state_pairs[0]
        return state
    total_weight = sum(weight for weight, state in weight_state_pairs)
    assert total_weight > 0
    r = random.uniform(0, total_weight)
    cumulative_weight = 0
    for weight, state in weight_state_pairs:
        if weight < 0:
            raise ValueError('weights must be non-negative')
        cumulative_weight += weight
        if r < cumulative_weight:
            return state
    raise ValueError('no choice was made')

def select(data, n):
    """
    This is from ActiveState recipe 269554.
    This is not by me.
    Find the nth rank ordered element (the least value has rank 0).
    """
    data = list(data)
    if not 0 <= n < len(data):
        raise ValueError('not enough elements for the given rank')
    while True:
        pivot = random.choice(data)
        pcount = 0
        under, over = [], []
        uappend, oappend = under.append, over.append
        for elem in data:
            if elem < pivot:
                uappend(elem)
            elif elem > pivot:
                oappend(elem)
            else:
                pcount += 1
        if n < len(under):
            data = under
        elif n < len(under) + pcount:
            return pivot
        else:
            data = over
            n -= len(under) + pcount

class Cache:
    def __init__(self, callback, cache_limit):
        self.callback = callback
        self.cache_limit = cache_limit
        self.cache = {}
    def __call__(self, arg):
        try:
            return self.cache[arg]
        except KeyError:
            value = self.callback(arg)
            if not self.is_full():
                self.cache[arg] = value
            return value
    def is_full(self):
        if self.cache_limit is None:
            return False
        return (len(self.cache) >= self.cache_limit)


class TestUtil(unittest.TestCase):

    def test_chopped_bresenham(self):
        """
        This example was inspired by wikipedia.
        U{image<http://en.wikipedia.org/wiki/Image:Bresenham.svg>}
        """
        seq = (1, 1, 1, 2, 2, 3, 0, 0, 0, 0, 0)
        self.assertEquals(len(seq), 11)
        categories = tuple(bresenham_line(len(seq), 5))
        self.assertEquals(categories, (0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4))
        result = tuple(chopped_bresenham(seq, 5))
        self.assertEquals(result, ((1,1), (1,2), (2,3,0), (0,0), (0,0)))

    def test_select(self):
        unsorted_arr = (1, 5, 2, 7, 8, 9, 1)
        self.assertEquals(select(unsorted_arr, 3), 5)


if __name__ == '__main__':
    unittest.main()
