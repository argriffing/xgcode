"""Find centers of clusters of some binary strings.
"""

from StringIO import StringIO
import time
import itertools
import optparse
import random

from SnippetUtil import HandlingError
import Util
import Form
import Progress

g_matrix = [
        [0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,1,1,0,1,0,0,0],
        [0,0,1,0,1,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,1,1,0,1,1,1,0,0,1,1,1,0,0,1,0,1,1,1],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1],
        [0,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,1,0,0,0,0,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [1,1,0,1,0,0,1,0,1,0,1,1,0,0,0,1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,0,1],
        [0,0,1,0,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,0,0,1,1,0,1,0,1,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,1,1,0,1,0,0,1,0,1,0,0,0],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],
        [0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0],
        [0,0,1,0,1,1,1,1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,1,1,1],
        [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0],
        [0,0,1,1,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1],
        [0,0,1,0,1,1,1,1,1,0,0,0,1,1,0,0,0,1,0,1,0,1,0,1,1,1,0,0,0,1,0,0,0,1,0,1,1,1],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0],
        [0,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
        [1,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,0,0],
        [0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,0],
        [0,0,1,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1],
        [0,0,1,1,1,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,1,0,1,1,0]]

def get_hamming_distance(a, b):
    """
    @param a: a sequence of ones and zeros
    @param b: a sequence of ones and zeros
    @return: the hamming distance
    """
    return Util.hamming_distance(a, b)

def get_folded_distance(a, b):
    """
    @param a: a sequence of ones and zeros
    @param b: a sequence of ones and zeros
    @return: a distance
    """
    max_distance = len(a)
    hdist = Util.hamming_distance(a, b)
    return min(hdist, max_distance - hdist)

def get_form():
    """
    @return: the body of a form
    """
    bstrings = '\n'.join(''.join(str(x) for x in row) for row in g_matrix)
    form_objects = [
            Form.Integer('ncenters', 'find this many centers', 3, low=1),
            Form.RadioGroup('options', 'use this distance', [
                Form.RadioItem('folded', 'folded hamming distance', True),
                Form.RadioItem('hamming', 'hamming distance')]),
            Form.MultiLine('bstrings', 'binary strings', bstrings)]
    return form_objects

def get_centers(arrs, ncenters, metric, deadline):
    """
    @param arrs: binary arrays
    @param ncenters: the number of centers to find
    @param metric: this function measures the distance between points
    @param deadline: raise an error if it gets later than this
    @return: the maximum distance to a center, and the list of centers
    """
    assert ncenters > 0
    assert ncenters <= len(arrs)
    # create the distance matrix
    n = len(arrs)
    D = [[metric(a,b) for b in arrs] for a in arrs]
    # get the centers
    best_max_distance = None
    best_centers = None
    for center_indices in itertools.combinations(range(n), ncenters):
        if deadline and time.time() > deadline:
            raise HandlingError('im sorry this calculation has exceeded my attention span')
        max_distance = max(min(D[c][t] for c in center_indices) for t in range(n))
        if not best_centers or max_distance < best_max_distance:
            best_max_distance = max_distance
            best_centers = [arrs[i] for i in center_indices]
    return best_max_distance, best_centers

def get_result_string(arrs, max_distance, centers):
    """
    @param arrs: binary arrays
    @param max_distance: the max distance of a binary string to the nearest center
    @param centers: the centers found by the search
    """
    out = StringIO()
    print >> out, 'number of binary strings:', len(arrs)
    print >> out, 'length of each binary string:', len(arrs[0])
    print >> out, 'maximum distance from a center:', max_distance
    print >> out, 'centers:'
    for center in centers:
        print >> out, ''.join(str(x) for x in center)
    return out.getvalue().strip()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the bstring lines
    lines = StringIO(fs.bstrings).readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    # do some validation
    if not lines:
        raise HandlingError('expected some binary strings')
    for line in lines:
        if not set(line) <= set('01'):
            raise HandlingError('expected only binary strings')
    d = len(lines[0])
    for line in lines:
        if len(line) != d:
            raise HandlingError('expected each binary string to have the same length')
    # convert the binary strings to arrays of integers
    arrs = [[int(x) for x in line] for line in lines]
    # define the distance metric
    if fs.hamming:
        metric = get_hamming_distance
    elif fs.folded:
        metric = get_folded_distance
    # do the analysis
    deadline = time.time() + 5.0
    max_distance, centers = get_centers(arrs, fs.ncenters, metric, deadline)
    result_string = get_result_string(arrs, max_distance, centers)
    return [('Content-Type', 'text/plain')], result_string

def get_histogram(arr):
    """
    @param arr: a sequence of numbers
    """
    d = {}
    for value in arr:
        count = d.get(value, 0)
        d[value] = count + 1
    return d

if __name__ == '__main__':
    # parse the options
    parser = optparse.OptionParser()
    parser.add_option('--ncenters', dest='ncenters', default=3, type='int', help='number of centers')
    parser.add_option('--do-permutations', action='store_true', dest='do_permutations', default=False, help='do a permutation test')
    options, args = parser.parse_args()
    # set the options
    arrs = g_matrix
    ncenters = options.ncenters
    metric = get_folded_distance
    deadline = None
    # show a preliminary data summary
    narrs = len(arrs)
    print 'number of binary strings:', narrs
    print 'length of each binary string:', len(arrs[0])
    print
    # get the list of numbers of centers to evaluate
    ncenter_list = range(1, ncenters+1) + range(narrs-ncenters+1, narrs+1)
    # get results for the observed data
    for i in ncenter_list:
        max_distance, centers = get_centers(arrs, i, metric, deadline)
        print 'results for', i, ('centers:' if i > 1 else 'center:')
        print 'maximum distance from a center:', max_distance
        print 'centers:'
        for center in centers:
            print ''.join(str(x) for x in center)
        print
    # do some permutation tests
    if options.do_permutations:
        print 'doing some permutation tests'
        results = []
        n = 100
        pbar = Progress.Bar(n)
        for i in range(n):
            # shuffle the binary strings
            for arr in arrs:
                random.shuffle(arr)
            # get results for the shuffled strings
            row = []
            for j in ncenter_list:
                max_distance, centers = get_centers(arrs, j, metric, deadline)
                row.append(max_distance)
            results.append(row)
            pbar.update(i+1)
        pbar.finish()
        columns = zip(*results)
        print
        print 'permutation test results:'
        for j, column in zip(ncenter_list, columns):
            print 'max distance histogram for', j, ('centers:' if j > 1 else 'center:')
            d = get_histogram(column)
            for value, count in sorted(d.items()):
                print '%d:\t%d' % (value, count)
            print

