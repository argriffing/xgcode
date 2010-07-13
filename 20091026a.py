"""Examine the Brownian motion approximation of allele frequencies.

Examine topological properties of the Brownian motion approximation
of allele frequencies.
In "Phylogenetic Analysis: Models and Estimation Procedures",
Cavalli-Sforza and Edwards use transformation of allele frequencies
due to Fisher.
In their paper they look at a projection from a Euclidean space
that includes a time dimension
onto the space of the current time: the "now" plane"
I want to consider a more general projection
whose lost dimensions are not time,
but a separate dimension for each ancestral node.
To make the comparison more direct,
I will use allele frequencies as the data that generate this space.
"""


from StringIO import StringIO
import math
import random

import numpy as np

from SnippetUtil import HandlingError
import SnippetUtil
import Form
import FormOut
import NewickIO
import Euclid
import NeighborJoining
import FelTree

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def do_analysis():
    # initialize estimates of allele frequency angles among 4 populations for 5 alleles
    ZZ = (.0000, .0000, .0000, .0000, .0000)
    SB = (.2585, .7274, .2222, .5091, .0000)
    SE = (.1822, .4718, .2852, .2170, .0000)
    SK = (.1874, .1461, .3050, .2883, .1132)
    BE = (.1094, .4907, .1554, .2922, .0000)
    BK = (.2270, .6398, .1024, .7974, .1132)
    EK = (.2178, .4465, .2293, .5052, .1132)
    V_list = [
            [ZZ, SB, SE, SK],
            [SB, ZZ, BE, BK],
            [SE, BE, ZZ, EK],
            [SK, BK, EK, ZZ]]
    V = np.array([[np.array(x).var() for x in row] for row in V_list])
    return str(V)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response
    out = StringIO()
    print >> out, do_analysis()
    # write the response
    return [('Content-Type', 'text/plain')], out.getvalue().strip()

if __name__ == '__main__':
    print do_analysis()
