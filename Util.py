
import random
import unittest
import itertools
from math import log

import scipy.stats
from scipy.special import gammaln

def binomial_log_pmf(observed_n, max_n, p_success):
    #TODO special cases
    accum = 0
    accum += gammaln(max_n + 1)
    accum -= gammaln(observed_n + 1)
    accum -= gammaln((max_n - observed_n) + 1)
    accum += observed_n * log(p_success)
    accum += (max_n - observed_n) * log(1.0 - p_success)
    return accum

def geometric_log_pmf(observed_n, pr):
    """
    @param observed_n: the number of completed events
    @param pr: the probability of quitting
    """
    if pr == 0.0:
        return float('-inf')
    if pr == 1.0:
        if observed_n:
            return float('-inf')
        else:
            return log(pr)
    return observed_n * log(1.0 - pr) + log(pr)

def poisson_log_pmf(observed_n, expected_n):
    if not expected_n:
        if observed_n:
            return float('-inf')
        else:
            return 0.0
    accum = 0
    accum += observed_n * log(expected_n)
    accum -= expected_n
    accum -= gammaln(observed_n+1)
    return accum

def multinomial_log_pmf(distribution, counts):
    """
    This should be in scipy.stats but it isn't.
    @param distribution: the distribution over classes
    @param counts: the observed counts over classes
    """
    # check for a degeneracy
    for d, c in zip(distribution, counts):
        if c and not d:
            return float('-inf')
    n = sum(counts)
    # initialize the log probability mass
    accum = 0
    # add the contribution of n to the multinomial coefficient
    if n > 1:
        accum += gammaln(n+1)
    # add the contribution of the counts to the multinomial coefficient
    accum -= sum(gammaln(count+1) for count in counts if count > 1)
    # add the contribution of probabilities
    for p, count in zip(distribution, counts):
        if count:
            accum += count * log(p)
    return accum

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

def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
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

def stripped_lines(lines):
    """
    This function yields nonempty stripped lines.
    """
    for line in lines:
        line = line.strip()
        if line:
            yield line

def dot_product(va, vb):
    assert len(va) == len(vb)
    return sum(a * b for a, b in zip(va, vb))

def product(numbers):
    x = 1
    for number in numbers:
        x *= number
    return x

def rle(sequence):
    """
    Yield (value, length) pairs.
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
    Yield chunks of a sequence whose lengths differ from each other by at most one.
    @param sequence: can be iterated and its length can be found
    @param nchunks: the number of yielded chunks
    """
    assert len(sequence) >= nchunks
    chunk = []
    current_category = None
    for value, category in zip(sequence, bresenham_line(len(sequence), nchunks)):
        if category == current_category:
            chunk.append(value)
        else:
            if current_category is not None:
                if type(sequence) is str:
                    yield ''.join(chunk)
                else:
                    yield tuple(chunk)
            current_category = category
            chunk = [value]
    if chunk:
        if type(sequence) is str:
            yield ''.join(chunk)
        else:
            yield tuple(chunk)

def chopped(sequence, size):
    """
    Yields regular sized chunks of a sequence, but the last one may be ragged.
    Notice that strings are treated differently than other iterables.
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
    This is from ActiveState recipe 269554 (contrib).
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
            if (self.cache_limit is None) or (len(self.cache) < self.cache_limit):
                self.cache[arg] = value
            return value


class TestUtil(unittest.TestCase):

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

    def test_chopped_bresenham(self):
        """
        Inspired by this U{image<http://en.wikipedia.org/wiki/Image:Bresenham.svg>}.
        """
        seq = (1, 1, 1, 2, 2, 3, 0, 0, 0, 0, 0)
        self.assertEquals(len(seq), 11)
        categories = tuple(bresenham_line(len(seq), 5))
        self.assertEquals(categories, (0, 0, 1, 1, 2, 2, 2, 3, 3, 4, 4))
        result = tuple(chopped_bresenham(seq, 5))
        self.assertEquals(result, ((1,1), (1,2), (2,3,0), (0,0), (0,0)))

    def test_dot_product(self):
        self.assertEquals(dot_product((1, 2, 3), (4, 5, 6)), 1*4 + 2*5 + 3*6)

    def test_select(self):
        unsorted_arr = (1, 5, 2, 7, 8, 9, 1)
        self.assertEquals(select(unsorted_arr, 3), 5)

    def test_poisson_log_pmf(self):
        observed_n = 60
        expected_n = 20
        likelihood = scipy.stats.poisson.pmf(observed_n, expected_n)
        expected = log(likelihood)
        observed = poisson_log_pmf(observed_n, expected_n)
        self.assertAlmostEqual(expected, observed)

    def test_geometric_log_pmf(self):
        obs = 5
        pr = 0.1
        scipy_result = log(scipy.stats.geom.pmf(obs, pr, loc=-1))
        util_result = geometric_log_pmf(obs, pr)
        self.assertAlmostEqual(scipy_result, util_result)


if __name__ == '__main__':
    unittest.main()
