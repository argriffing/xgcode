r"""
Plot some scaled fixation probabilities assuming recessivity/dominance.

This script plots some values of the limit
\[
\lim_{N \to \infty}
2N \cdot u \left( p=\frac{1}{2N}, c=Ns, D=z \cdot \text{sgn} (s) \right)
\]
for various combinations of \(c\) and \(z\).
Most of that notation follows Kimura 1957.
When \(z=0\) this reduces to a classic approximation.
"""

from StringIO import StringIO
import math

import numpy

import Form
import FormOut
import kimrecessive
import MatrixUtil
import RUtil
from RUtil import mk_call_str

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.ImageFormat(),
            ]
    return form_objects

def get_form_out():
    return FormOut.Image('scaled-fixation-plot')

def get_response_content(fs):
    # create the R table string and scripts
    headers = [
            'z',
            'c.neg2.0',
            'c.neg0.5',
            'c.0.5',
            'c.2.0',
            #'c.a',
            #'c.b',
            #'c.c',
            #'c.d',
            ]
    #C = numpy.array([-0.5, -0.2, 0.2, 0.5], dtype=float)
    #C = numpy.array([-1.0, -0.4, 0.4, 1.0], dtype=float)
    C = numpy.array([-2.0, -0.5, 0.5, 2.0], dtype=float)
    Z = numpy.linspace(-3, 3, 101)
    # get the data for the R table
    arr = []
    for z in Z:
        row = [z]
        for c in C:
            rate = 1.0 / kimrecessive.denom_piecewise(c, z*numpy.sign(c))
            row.append(rate)
        arr.append(row)
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
    print >> out, mk_call_str('melt', 'my.table', id='"z"')
    print >> out, 'ggplot(data=my.table.long,'
    print >> out, mk_call_str('aes', x='z', y='value', colour='variable')
    print >> out, ') + geom_line()',
    print >> out, '+',
    print >> out, mk_call_str(
            'xlim',
            mk_call_str('min', 'my.table.long$z'),
            mk_call_str('max', 'my.table.long$z'),
            ),
    print >> out, '+',
    print >> out, mk_call_str(
            'ylim', '0',
            mk_call_str('max', 'my.table.long$value'))
    return out.getvalue()

