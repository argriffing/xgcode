"""
Scatter plot Shannon entropy vs. an entropy-like function.

Shannon entropy is negative expectation of prob log prob.
An analog is log ratio sum of squares to sum of cubes.
"""

from StringIO import StringIO
import random
import math

import numpy as np

import Form
import FormOut
import RUtil
import mrate
import ctmcmi

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('mutual-info-and-analog')

def get_response_content(fs):
    # create the R table string and scripts
    headers = [
            'entropy',
            'analog']
    distributions = []
    nstates = 4
    npoints = 5000
    arr = []
    best_pair = None
    for i in range(npoints):
        weights = [random.expovariate(1) for j in range(nstates)]
        total = sum(weights)
        distn = [x / total for x in weights]
        entropy = -sum(p * math.log(p) for p in distn)
        sum_squares = sum(p*p for p in distn)
        sum_cubes = sum(p*p*p for p in distn)
        analog = math.log(sum_squares / sum_cubes)
        row = [entropy, analog]
        arr.append(row)
        dist = (entropy - 1.0)**2 + (analog - 0.4)**2
        if (best_pair is None) or (dist < best_pair[0]):
            best_pair = (dist, distn)
    # get the R table
    table_string = RUtil.get_table_string(arr, headers)
    # get the R script
    out = StringIO()
    title = ', '.join(str(x) for x in best_pair[1])
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$entropy',
            'my.table$analog',
            pch='20',
            main='"%s"' % title)
    script = out.getvalue()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

