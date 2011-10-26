"""Look at MDS projections of two inter-city distances.
"""

from StringIO import StringIO
import os
import math

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Util
import Form
import FormOut
import GPS
import Euclid
import const

g_mds_data = const.read('20100119b')

def parse_lines(lines):
    """
    The input lines have a special format.
    The first nonempty line is a header.
    The subsequent lines are whitespace separated values.
    The first value is the city name.
    The next four values are latitude and longitude minutes and degrees.
    @param lines: stripped input lines
    @return: (city, lat_deg, lat_min, lon_deg, lon_min) tuples
    """
    lines = [line for line in lines if line]
    if not lines:
        raise HandlingError('no input was found')
    if len(lines) < 2:
        raise HandlingError('expected at least one header and data line')
    result = []
    for line in lines[1:]:
        values = line.split()
        if len(values) != 5:
            raise HandlingError('expected five values per data line')
        city, latd, latm, lond, lonm = values
        try:
            latd = float(latd)
            latm = float(latm)
            lond = float(lond)
            lonm = float(lonm)
        except ValueError, e:
            raise HandlingError('error reading a value as a number')
        row = (city, latd, latm, lond, lonm)
        result.append(row)
    return result

def get_form():
    """
    @return: a list of form objects
    """
    # define the list of form objects
    lines = [line.strip() for line in g_mds_data.splitlines()]
    form_objects = [
            Form.MultiLine('datalines', 'locations', '\n'.join(lines))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the lat-lon points from the input
    lines = Util.get_stripped_lines(fs.datalines.splitlines())
    rows = parse_lines(lines)
    latlon_points = []
    city_names = []
    for city, latd, latm, lond, lonm in rows:
        lat = math.radians(GPS.degrees_minutes_to_degrees(latd, latm))
        lon = math.radians(GPS.degrees_minutes_to_degrees(lond, lonm))
        latlon_points.append((lat, lon))
        city_names.append(city)
    npoints = len(latlon_points)
    # start writing the response
    np.set_printoptions(linewidth=200)
    out = StringIO()
    radius = GPS.g_earth_radius_miles
    for dfunc, name in (
            (GPS.get_arc_distance, 'great arc'),
            (GPS.get_euclidean_distance, 'euclidean')):
        # define the edm whose elements are squared euclidean-like distances
        edm = np.zeros((npoints, npoints))
        D = np.zeros((npoints, npoints))
        for i, pointa in enumerate(latlon_points):
            for j, pointb in enumerate(latlon_points):
                D[i, j] = dfunc(pointa, pointb, radius)
                edm[i, j] = D[i, j]**2
        print >> out, name, 'distances:'
        print >> out, D
        print >> out
        print >> out, name, 'EDM:'
        print >> out, edm
        print >> out
        G = Euclid.edm_to_dccov(edm)
        print >> out, name, 'Gower centered matrix:'
        print >> out, G
        print >> out
        spectrum = np.array(list(reversed(sorted(np.linalg.eigvals(G)))))
        print >> out, name, 'spectrum of Gower centered matrix:'
        for x in spectrum:
            print >> out, x
        print >> out
        print >> out, name, 'rounded spectrum:'
        for x in spectrum:
            print >> out, '%.1f' % x
        print >> out
        mds_points = Euclid.edm_to_points(edm)
        print >> out, '2D MDS coordinates:'
        for name, mds_point in zip(city_names, mds_points):
            x = mds_point[0]
            y = mds_point[1]
            print >> out, '\t'.join(str(x) for x in [name, x, y])
        print >> out
        # break between distance methods
        print >> out
    # return the response
    return out.getvalue()
