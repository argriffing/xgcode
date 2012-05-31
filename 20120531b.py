r"""
Plot max mutual information for several selection models and their limits.

The mutation process is site-independent 3-site 2-state-per-site.
"""

from StringIO import StringIO
import argparse
import math
import time
import random
import heapq
from itertools import product

import numpy as np
import scipy
from scipy import linalg
from scipy import optimize

import Form
import FormOut
import ctmcmi
import mrate
import divtime
import cheeger
import MatrixUtil
from MatrixUtil import ndot
import RUtil
from RUtil import mk_call_str
import evozoo

# variable name, description, python object
g_process_triples = [
        ('cube', '3d cube',
            evozoo.Hypercube_2_3_0()),
        ('cycle', 'maximal induced cycle',
            evozoo.Coil_2_3_0()),
        ('path', 'maximal induced path',
            evozoo.Snake_2_3_0()),
        ]

def get_form():
    """
    @return: the body of a form
    """
    check_items = [Form.CheckItem(a, b, True) for a, b, c in g_process_triples]
    return [
            Form.CheckGroup('logs', 'plot process reducibility limits', [
                Form.CheckItem(
                    'log4', 'log(4) cube limit', True),
                Form.CheckItem(
                    'log3', 'log(3) cycle and path limit', True)]),
            Form.CheckGroup(
                'processes', 'plot mutual information', check_items),
            Form.Float('start_time', 'start time', '0.04', low_exclusive=0),
            Form.Float('stop_time', 'stop time', '0.1', low_exclusive=0),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def gen_overdispersed_events(low, high):
    """
    This is a generator that samples overdispersed events in an interval.
    It is useful for plotting.
    """
    # Create the min heap.
    # The triples are the neg length, the low, and the high values.
    # Popping from the queue will return the longest interval.
    # The queue grows linearly with the number of sampled events.
    yield low
    yield high
    q = [(low-high, low, high)]
    while True:
        dummy, a, b = heapq.heappop(q)
        mid = random.uniform(a, b)
        heapq.heappush(q, (a-mid, a, mid))
        heapq.heappush(q, (mid-b, mid, b))
        yield mid

class OptDep:
    def __init__(self, zoo_obj, t, f_info):
        """
        @param zoo_obj: an object from the evozoo module
        @param t: divergence time
        @param f_info: info function that takes a rate matrix and a time
        """
        self.zoo_obj = zoo_obj
        self.t = t
        self.f_info = f_info
    def __call__(self, X):
        """
        @param X: some log ratio probabilities
        @return: neg info value for minimization
        """
        distn = self.zoo_obj.get_distn(X)
        Q = self.zoo_obj.get_rate_matrix(X)
        return -self.f_info(Q, distn, self.t)

def get_response_content(fs):
    # hardcode the amount of time allowed
    nseconds = 8
    # validate and store user input
    if fs.stop_time <= fs.start_time:
        raise ValueError('check the start and stop times')
    f_info = ctmcmi.get_mutual_info_known_distn
    requested_triples = []
    for triple in g_process_triples:
        name, desc, zoo_obj = triple
        if getattr(fs, name):
            requested_triples.append(triple)
    # define the R table headers
    headers = ['t']
    if fs.log4:
        headers.append('log.4')
    if fs.log3:
        headers.append('log.3')
    r_names = [a.replace('_', '.') for a, b, c in requested_triples]
    headers.extend(r_names)
    # Spend a lot of time doing the optimizations
    # to construct the points for the R table.
    ntimes = 100
    incr = (fs.stop_time - fs.start_time) / float(ntimes - 1)
    times = [fs.start_time + i*incr for i in range(ntimes)]
    arr = []
    for t in times:
        row = [t]
        if fs.log4:
            row.append(math.log(4))
        if fs.log3:
            row.append(math.log(3))
        for python_name, desc, zoo_obj in requested_triples:
            X = np.array([])
            info_value = f_info(
                    zoo_obj.get_rate_matrix(X),
                    zoo_obj.get_distn(X),
                    t)
            row.append(info_value)
        arr.append(row)
    # create the R table string and scripts
    # get the R table
    table_string = RUtil.get_table_string(arr, headers)
    # get the R script
    script = get_ggplot()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

def get_ggplot():
    out = StringIO()
    print >> out, mk_call_str('require', '"reshape"')
    print >> out, mk_call_str('require', '"ggplot2"')
    print >> out, 'my.table.long <-',
    print >> out, mk_call_str('melt', 'my.table', id='"t"')
    print >> out, 'ggplot(data=my.table.long,'
    print >> out, mk_call_str('aes', x='t', y='value', colour='variable')
    print >> out, ') + geom_line()',
    print >> out, '+',
    print >> out, mk_call_str(
            'xlim',
            mk_call_str('min', 'my.table.long$t'),
            mk_call_str('max', 'my.table.long$t')),
    print >> out, '+',
    print >> out, mk_call_str(
            'ylim', '0',
            mk_call_str('max', 'my.table.long$value'))
    return out.getvalue()

