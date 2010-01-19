"""
This module is about latitude and longitude.

Haversine formula and unit test is from a snippet by Wayne Dyck.
Cartesian coordinate conversion and the associated unit test
are from a gps webpage at Frank Wattenberg's site
called Spherical Coordinates and the GPS.
"""

import unittest
import math

import numpy as np


g_const_data = 'const_data'

g_miles_per_kilometer = 0.621371192237334
g_kilometers_per_mile = 1 / g_miles_per_kilometer

g_earth_radius_kilometers = 6371.0
g_earth_radius_miles = g_earth_radius_kilometers * g_miles_per_kilometer


def minutes_seconds_to_minutes(m, s):
    """
    @param m: minutes
    @param s: seconds
    @return: minutes
    """
    return float(m) + float(s) / 60

def degrees_minutes_to_degrees(d, m):
    """
    @param d: degrees
    @param m: minutes
    @return: degrees
    """
    return float(d) + float(m) / 60

def get_arc_distance(origin_latlon, target_latlon, radius):
    """
    Return the great arc distance between two points.
    @param origin: (latitude_north, longitude_west) radians of the origin
    @param target: (latitude_north, longitude_west) radians of the target
    @param radius: the radius of the sphere
    @return: great arc distance in the units of the radius of the sphere
    """
    src_lat, src_lon = origin_latlon
    dst_lat, dst_lon = target_latlon
    dlat = dst_lat - src_lat
    dlon = dst_lon - src_lon
    sin_dlat_d2_2 = math.sin(dlat/2)**2
    sin_dlon_d2_2 = math.sin(dlon/2)**2
    a = sin_dlat_d2_2 + math.cos(src_lat) * math.cos(dst_lat) * sin_dlon_d2_2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1.0 - a))
    d = radius * c
    return d

def get_euclidean_distance(origin_latlon, target_latlon, radius):
    """
    Return the Euclidean distance between two (lat, lon) points.
    The returned distance is the straight-line-through-the-earth distance.
    @param origin: (latitude_north, longitude_west) radians of the origin
    @param target: (latitude_north, longitude_west) radians of the target
    @param radius: the radius of the sphere
    @return: great arc distance in the units of the radius of the sphere
    """
    origin = latlon_to_euclidean(origin_latlon, radius)
    target = latlon_to_euclidean(target_latlon, radius)
    return np.linalg.norm(target - origin)

def latlon_to_euclidean(point, radius):
    """
    @param point: (latitude_north, longitude_west) radians of the point
    @param radius: the radius of the sphere
    @return: point in 3d euclidean space
    """
    lat, lon = point
    phi = (math.pi / 2) - lat
    theta = -lon
    x = radius * math.cos(theta) * math.sin(phi)
    y = radius * math.sin(theta) * math.sin(phi)
    z = radius * math.cos(phi)
    return np.array([x, y, z])


class TestGPS(unittest.TestCase):

    def test_earth_radius(self):
        self.assertTrue(3500 < g_earth_radius_miles < 4500)
        self.assertTrue(5500 < g_earth_radius_kilometers < 7000)

    def test_get_arc_distance(self):
        seattle_lat = math.radians(47.621800)
        seattle_lon = math.radians(-122.350326)
        seattle = (seattle_lat, seattle_lon)
        olympia_lat = math.radians(47.041917)
        olympia_lon = math.radians(-122.893766)
        olympia = (olympia_lat, olympia_lon)
        radius = g_earth_radius_kilometers
        expected = 76.386615799548693
        observed = get_arc_distance(seattle, olympia, radius)
        self.assertAlmostEqual(expected, observed)

    def test_latlon_to_euclidean(self):
        """
        Test the conversion from (lat, lon) to (x, y, z).
        This function uses a custom constant for the Earth's radius.
        """
        lat = math.radians(degrees_minutes_to_degrees(35, 17.299))
        lon = math.radians(degrees_minutes_to_degrees(120, 39.174))
        observed = latlon_to_euclidean((lat, lon), 6367).round()
        expected = np.array([-2650, -4471, 3678])
        self.assertTrue(np.allclose(observed, expected))

    def test_foo(self):
        seattle_lat = math.radians(47.621800)
        seattle_lon = math.radians(-122.350326)
        seattle = (seattle_lat, seattle_lon)
        olympia_lat = math.radians(47.041917)
        olympia_lon = math.radians(-122.893766)
        olympia = (olympia_lat, olympia_lon)
        radius = g_earth_radius_kilometers
        print get_arc_distance(seattle, olympia, radius)
        print get_euclidean_distance(seattle, olympia, radius)


if __name__ == '__main__':
    unittest.main()
