"""
Estimate mutual intersection times among piecewise Bezier paths.

The idea is that we start with a bunch of piecewise Bezier curves
and we want to find their mutual intersection times,
but we do not care about either their self intersection times
or the identities of the curves which they intersect.
Also we might only care about their intersections under some transformation.
"""

from StringIO import StringIO
from collections import defaultdict

import numpy as np

import bezier
import pcurve


class _OwnedBezierChunk(bezier.BezierChunk):
    """
    Allow Bezier chunks to recognize member of the same parent curve.
    This allows recognition even when the chunks are mixed
    together in an amorphous soup for mutual intersection detection.
    """
    def __init__(self, *args, **kwargs):
        bezier.BezierChunk.__init__(self, *args, **kwargs)
        self.parent_ref = None
    def split(self, t):
        a, b = bezier.BezierChunk.split(self, t)
        a.parent_ref = self.parent_ref
        b.parent_ref = self.parent_ref
        return a, b


def get_intersection_times(bpaths, trans, min_gridsize, min_spatial_gap):
    """
    Return lists of filtered intersection times.
    Points are represented by numpy arrays.
    @param bpaths: a sequence of pcurve.BezierPath objects
    @param trans: a point to point transformation
    @param min_gridsize: a grid resolution for intersection detection error
    @param min_spatial_gap: a resolution for spacing between intersections
    @return: a sequence, conformant with bpaths, of lists of intersection times
    """
    # Clone the bpaths in preparation for transformation.
    # It is not enough to create transformed cloned bchunks
    # because the bpaths must be used for time filtering.
    transformed_bpaths = []
    for bpath in bpaths:
        transformed_bpath = bpath.clone()
        transformed_bpath.transform(trans)
        transformed_bpaths.append(transformed_bpath)
    # create a soup of owned transformed bchunks
    bchunk_soup = []
    for i, bpath in enumerate(transformed_bpaths):
        for b in bpath.bchunks:
            bnew = _OwnedBezierChunk(
                    b.start_time, b.stop_time, b.p0, b.p1, b.p2, b.p3)
            bnew.parent_ref = i
            bchunk_soup.append(bnew)
    # get the list of refined bchunks that intersect with something interesting
    intersecting_bchunks = _get_bchunk_intersections(bchunk_soup, min_gridsize)
    # get the set of times per bpath index
    unfiltered_time_map = defaultdict(set)
    for b in intersecting_bchunks:
        unfiltered_time_map[b.parent_ref].update((b.start_time, b.stop_time))
    # filter the times
    filtered_times_list = []
    for i, bpath in enumerate(transformed_bpaths):
        unfiltered_times = unfiltered_time_map.get(i, [])
        filtered_times = _filter_intersection_times(
                bpath, unfiltered_times, min_spatial_gap)
        filtered_times_list.append(filtered_times)
    return filtered_times_list

def _filter_intersection_times(bpath, raw_intersection_times, min_spatial_gap):
    """
    Collapse intersection clusters.
    @param bpath: a pcurve.BezierPath object
    @param raw_intersection_times: a collection of intersection times
    @param min_spatial_gap: minimium spatial gap between sequential events
    @return: filtered time sequence
    """
    # first sort the intersection times
    times = sorted(raw_intersection_times)
    # group together times that are indistinguishable
    groups = []
    last_point = None
    g = []
    for t in times:
        point = bpath.evaluate(t)
        # if we have seen a previous point then check the gap
        if g:
            gap = np.linalg.norm(point - last_point)
            # if the gap is large then start a new group
            if gap >= min_spatial_gap:
                groups.append(g)
                g = []
        # append the current time to the current group
        g.append(t)
        # remember the most recent point
        last_point = point
    if g:
        groups.append(g)
    # return the sequence of group midpoints
    return [0.5 * (g[0] + g[-1]) for g in groups]

def _get_bchunk_intersections(bchunks, min_gridsize):
    """
    This is essentially a dumb search.
    It looks for collisions of smaller and smaller curve pieces
    on finer and finer grids.
    The only smartness is that if a large curve piece
    has no collisions on a coarse grid,
    then it is not subdivided for consideration in a finer grid search.
    Self intersections of curves are not considered.
    @param bchunks: a collection of OwnedBezierChunk objects
    @param min_gridsize: a float lower bound resolution
    @return: a collection of refined intersecting bchunks
    """
    too_many_things = 1e6
    # Maintain the invariant that potentially intersecting chunks
    # have a diameter of no more than twice the gridsize.
    gridsize = 0.5 * max(b.get_diameter() for b in bchunks)
    while len(bchunks) < too_many_things:
        # map each grid point to a set of nearby parent curves
        gridmap = defaultdict(set)
        for b in bchunks:
            for gridpoint in b.gen_bb_gridpoints(gridsize):
                gridmap[gridpoint].add(b.parent_ref)
        # Get the set of indices of bchunks
        # whose bounding boxes contain contested grid points.
        index_set = set()
        for i, b in enumerate(bchunks):
            for gridpoint in b.gen_bb_gridpoints(gridsize):
                if len(gridmap[gridpoint]) > 1:
                    index_set.add(i)
                    break
        # Cut the gridsize in half.
        gridsize *= 0.5
        # If the gridsize is below the min
        # then return the bchunks involved in putative intersections.
        if gridsize < min_gridsize:
            return [bchunks[index] for index in index_set]
        # Classify each intersecting bchunk by its diameter.
        bchunks_small = []
        bchunks_large = []
        for index in index_set:
            b = bchunks[index]
            if b.get_diameter() <= gridsize*2:
                bchunks_small.append(b)
            else:
                bchunks_large.append(b)
        # Iteratively bisect large intersecting chunks.
        while bchunks_large:
            b = bchunks_large.pop()
            for child in b.bisect():
                if child.get_diameter() <= gridsize*2:
                    bchunks_small.append(child)
                else:
                    bchunks_large.append(child)
        bchunks = bchunks_small
    out = StringIO()
    print >> out, 'too many things:'
    print >> out, len(bchunks), 'bchunks'
    print >> out, gridsize, 'gridsize'
    print >> out, min_gridsize, 'min_gridsize'
    raise ValueError(out.getvalue())

