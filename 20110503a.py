"""Check the variance of the difference from the leaf mean.

This is obsolete.
Do not use.
"""

from collections import defaultdict
from StringIO import StringIO
import random
import math

import Form
import FormOut
import Harmonic
import Ftree
import FtreeIO

g_default_tree = '((1:1, 2:0.5)6:1, (3:0.333333333333, 4:0.5)7:1, 5:1)8;'

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', g_default_tree)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def bridge(a, b, npillars, distance):
    """
    Sample a value at each pillar of a brownian bridge.
    The pillars are assumed to be equally spaced.
    @param a: Brownian value at the initial point
    @parma b: Brownian value at the final point
    @param npillars: the number of intervening points to sample
    @param distance: the distance betweeen the initial and final points
    @return: a list of samples
    """
    samples = []
    incr = distance / float(npillars + 1)
    prev = a
    for i in range(npillars):
        t = incr
        t2 = distance - incr*i
        mu = prev + (t / t2) * (b - prev)
        var = t*(t2 - t) / t2
        x = random.gauss(mu, math.sqrt(var))
        samples.append(x)
        prev = x
    return samples

def sample_brownian_motion(R, B):
    """
    Sample brownian motion on a tree.
    @param R: directed tree
    @param B: branch lengths
    @return: map from vertex to sample
    """
    r = Ftree.R_to_root(R)
    v_to_sample = {r : 0}
    v_to_sinks = Ftree.R_to_v_to_sinks(R)
    for v in Ftree.R_to_preorder(R):
        for sink in v_to_sinks[v]:
            u_edge = frozenset((v, sink))
            mu = v_to_sample[v]
            var = B[u_edge]
            v_to_sample[sink] = random.gauss(mu, math.sqrt(var))
    return v_to_sample

def get_samples(point_triples):
    """
    Return two dictionaries.
    The first dictionary maps a vertex to a vector of samples.
    The vertices in this first dictionary include
    all of the points of articulation in addition to all of the leaves.
    The second dictionary maps a (source, sink) vector pair
    defining a directed edge to a (distance, vector) pair.
    The distance in this second dictionary is the distance
    from the source vertex to the sampling point,
    and the vector is the vector of samples at that sampling point.
    @param foo:
    @return: vertex vectors and edge vectors
    """
    pass


def get_variances(T, B, N, nsamples, blen_granularity):
    """
    Look at the variance of the difference from the leaf mean.
    Do this simultaneously for multiple points
    on multiple branches on the tree.
    Report quadruples (var, source, sink, distance)
    where var is the estimated variance,
    source is the source vertex of an edge,
    sink is the sink vertex of an edge,
    and distance is the distance from the source.
    @param T: topology
    @param B: branch lengths
    @param nsamples: use this many samples per point
    @param blen_granularity: rough size of smallest segments
    @return: sequence of quadruples
    """
    pass

def get_response_content(fs):
    # read the tree
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    leaves = Ftree.T_to_leaves(T)
    internal = Ftree.T_to_internal_vertices(T)
    # root arbitrarily
    R = Ftree.T_to_R_canonical(T)
    # init some sampling parameters
    nsamples = 1000
    npillars = 10
    # Init the accumulators.
    # Accumulate the sum of squares of differences
    # and the sum of differences.
    # The differences are from the leaf mean.
    dsum = defaultdict(float)
    dsumsq = defaultdict(float)
    # Repeatedly sample using Brownian motion on the tree.
    for i in range(nsamples):
        # Sample using Brownian motion at vertices on the tree.
        v_to_sample = sample_brownian_motion(R, B)
        # Compute the mean at the leaves.
        mu = sum(v_to_sample[v] for v in leaves) / len(leaves)
        # Accumulate difference moments at vertices of the tree.
        for v, x in v_to_sample.items():
            dsum[(v, -1, -1)] += x-mu
            dsumsq[(v, -1, -1)] += (x-mu)**2
        # Sample using Brownian bridge on edges.
        for d_edge in R:
            u_edge = frozenset(d_edge)
            va, vb = d_edge
            a = v_to_sample[va]
            b = v_to_sample[vb]
            samples = bridge(a, b, npillars, B[u_edge])
            for i, x in enumerate(samples):
                dsum[(va, vb, i)] += x-mu
                dsumsq[(va, vb, i)] += (x-mu)**2
    quad = min((val, va, vb, i) for (va, vb, i), val in dsumsq.items())
    val, va, vb, i = quad
    # write the report
    out = StringIO()
    if i < 0:
        print >> out, 'min sum of squares was at vertex', N[va]
    else:
        print >> out, 'min sum of squares was at edge',
        print >> out, N[va], '--[', i, ']-->', N[vb]
    print >> out
    print >> out, 'the min sum of squares value was', val
    return out.getvalue()
