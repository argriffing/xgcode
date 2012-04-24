"""
Plot log(ratio(E(L))) and E(log(ratio(L))) for non-uniform stationary process.
"""

from StringIO import StringIO
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
            't',
            'mutual.info.mut',
            'mutual.info.analog.mut',
            'mutual.info.mutsel',
            'mutual.info.analog.mutsel']
    npoints = 100
    t_low = 0.001
    t_high = 5.0
    t_incr = (t_high - t_low) / (npoints - 1)
    t_values = [t_low + t_incr*i for i in range(npoints)]
    # get the mutation rate matrix
    """
    M = (1.0 / 2.0) * np.array([
        [0, 1, 0, 1],
        [1, 0, 1, 0],
        [0, 1, 0, 1],
        [1, 0, 1, 0]], dtype=float)
    """
    M = (1.0 / 3.0) * np.array([
        [0, 1, 1, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [1, 1, 1, 0]], dtype=float)
    M -= np.diag(np.sum(M, axis=1))
    # adjust the mutation process to give it new stationary probabilities
    #mut_weights = np.array([1, 1, 0.001, 0.001], dtype=float)
    mut_weights = np.array([1, 0.001, 1, 0.001], dtype=float)
    mut_distn = [x / np.sum(mut_weights) for x in mut_weights]
    M = mrate.to_gtr_halpern_bruno(M, mut_distn)
    # define the target stationary distribution
    p = math.log(2)
    q = 1 - p
    mutsel_weights = np.array([p, q/3, q/3, q/3], dtype=float)
    mutsel_distn = [x / np.sum(mutsel_weights) for x in mutsel_weights]
    # get the mutation selection balance rate matrix
    R = mrate.to_gtr_halpern_bruno(M, mutsel_distn)
    # get the data for the R table
    arr = []
    for t in t_values:
        mi_mut = ctmcmi.get_mutual_information(M, t)
        mi_analog_mut = ctmcmi.get_ll_ratio_wrong(M, t)
        mi_mutsel = ctmcmi.get_mutual_information(R, t)
        mi_analog_mutsel = ctmcmi.get_ll_ratio_wrong(R, t)
        row = [t, mi_mut, mi_analog_mut, mi_mutsel, mi_analog_mutsel]
        arr.append(row)
    # get the R table
    table_string = RUtil.get_table_string(arr, headers)
    # get the R script
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$mutual.info.mutsel',
            type='"n"',
            xlab='"time"',
            ylab='""',
            ylim=RUtil.mk_call_str('c', '0', '1.0'),
            main='"MI (black) and MI analog (red) over time"')
    print >> out, RUtil.mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.mut',
            col='"black"')
    print >> out, RUtil.mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.analog.mut',
            col='"red"')
    print >> out, RUtil.mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.mutsel',
            col='"black"')
    print >> out, RUtil.mk_call_str(
            'lines',
            'my.table$t',
            'my.table$mutual.info.analog.mutsel',
            col='"red"')
    script = out.getvalue()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

