"""Sample variances of differences of brownian motion on a graph.
"""

from StringIO import StringIO
import time
import random
import math
import optparse

from SnippetUtil import HandlingError
import Util
import MatrixUtil
import Form
import FormOut
import Progress

class Edge:
    def __init__(self, source, sink, distance):
        self.source = source
        self.sink = sink
        self.distance = distance
    def __str__(self):
        return ' '.join((self.source, self.sink, str(self.distance)))

g_default_edges = [
        Edge('a', 'b', 5.65),
        Edge('a', 'c', 6.64705882),
        Edge('a', 'd', 14.125),
        Edge('b', 'c', 22.6),
        Edge('b', 'd', 12.5555555555),
        Edge('c', 'd', 56.5)]

def get_form():
    """
    @return: the body of a form
    """
    default_edge_lines = [str(edge) for edge in g_default_edges]
    graphlines_data = '\n'.join(default_edge_lines)
    # define the form objects
    form_objects = [
            Form.MultiLine('graphlines',
                'edge sources, sinks, and distances', graphlines_data),
            Form.Float('epsilon',
                'use this much sloppiness in the conditioning', '1.0'),
            Form.Integer('nsamples',
                'accepted joint samples for variance estimation', '100')]
    return form_objects

def get_form_out():
    return FormOut.Report()

class RejectionError(Exception):
    pass

def get_joint_sample(edges, epsilon):
    """
    @param edges: a list of edge objects that define the graph
    @param epsilon: allow this much sloppiness in the conditioning
    @return: a dictionary mapping each state to a value
    """
    state_to_value = {}
    for edge in edges:
        if edge.source not in state_to_value:
            state_to_value[edge.source] = 0
        source_value = state_to_value[edge.source]
        putative_sink_value = random.gauss(source_value, math.sqrt(edge.distance))
        if edge.sink in state_to_value:
            error = abs(state_to_value[edge.sink] - putative_sink_value)
            if error > epsilon:
                raise RejectionError()
        else:
            state_to_value[edge.sink] = putative_sink_value
    return state_to_value

def get_variance(arr):
    """
    @param arr: an array of numbers
    @return: a variance
    """
    mu = sum(arr) / len(arr)
    return sum((value - mu)**2 for value in arr) / len(arr)

def get_zero_mean_variance(arr):
    """
    @param arr: an array of numbers
    @return: a variance
    """
    mu = 0
    return sum((value - mu)**2 for value in arr) / len(arr)

def process(edges, epsilon, nsamples, deadline, pbar):
    """
    @param edges: an ordered list of edge objects
    @param epsilon: this is how sloppy we can be
    @param nsamples: use this many accepted samples as the basis for the variance
    @param deadline: we raise an error this many seconds after 1970
    @param pbar: progress bar or None
    @return: text explaining our accomplishments
    """
    # get the set of states
    state_set = set()
    for edge in edges:
        state_set.add(edge.source)
        state_set.add(edge.sink)
    ordered_states = list(sorted(state_set))
    # initialize the pairwise distance lists
    pairwise_differences = {}
    for a in ordered_states:
        for b in ordered_states:
            pairwise_differences[(a,b)] = []
    # for each accepted joint sample get all pairwise distances
    naccepted_samples = 0
    while naccepted_samples < nsamples:
        if deadline and time.time() > deadline:
            raise HandlingError('sorry this was taking way too long')
        try:
            state_to_value = get_joint_sample(edges, epsilon)
        except RejectionError as e
            continue
        naccepted_samples += 1
        if pbar:
            pbar.update(naccepted_samples)
        for a in ordered_states:
            for b in ordered_states:
                difference = state_to_value[b] - state_to_value[a]
                pairwise_differences[(a,b)].append(difference)
    # write the variance of each pairwise distance
    out = StringIO()
    nstates = len(ordered_states)
    for i in range(nstates-1):
        for j in range(i+1, nstates):
            a = ordered_states[i]
            b = ordered_states[j]
            v = get_zero_mean_variance(pairwise_differences[(a,b)])
            print >> out, '\t'.join((a, b, str(v)))
    return out

def get_response_content(fs):
    # read the graph lines
    lines = Util.get_stripped_lines(StringIO(fs.graphlines))
    edges = []
    for line in lines:
        items = line.split()
        if len(items) != 3:
            raise HandlingError('expected three items per edge')
        source, sink, distance_string = items
        try:
            distance = float(distance_string)
        except ValueError as e
            raise HandlingError('expected each edge distance to be readable as a floating point number')
        edges.append(Edge(source, sink, distance))
    # if we exceed this many seconds since 1970 then we fail
    deadline = time.time() + 10
    # get the results
    out = process(edges, fs.epsilon, fs.nsamples, deadline, None)
    # return the response
    return out.getvalue()

def main(options):
    """
    @param options: parsed from the command line
    """
    edges = g_default_edges
    deadline = None
    pbar = Progress.Bar(options.nsamples)
    out = process(edges, options.epsilon, options.nsamples, deadline, pbar)
    pbar.finish()
    print out.getvalue().strip()

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('--nsamples', dest='nsamples', type='int', default=1000, help='number of samples')
    parser.add_option('--epsilon', dest='epsilon', type='float', default=0.5, help='sloppiness')
    options, args = parser.parse_args()
    main(options)

