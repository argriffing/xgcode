"""
Plot something that is like an information criterion from an email.

This is the plot of information criteria of two Jukes-Cantor
processes against each other.
"""

from StringIO import StringIO
import math

import Form
import FormOut
import RUtil

TOPLEFT = 'topleft'
BOTTOMLEFT = 'bottomleft'
TOPRIGHT = 'topright'
BOTTOMRIGHT = 'bottomright'

# wolframalpha plot:
# plot log(1 + 3*exp(-8/3 * x)), log(1 + exp(-4/3 * x)) , x=0 to 5

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nstates_a', 'first process state space size',
                4, low=2, high=9),
            Form.Integer('nstates_b', 'second process state space size',
                2, low=2, high=9),
            Form.Float('mu_a', 'first process randomization rate',
                '1.333333333', low_exclusive=0),
            Form.Float('mu_b', 'second process randomization rate',
                '0.666666666', low_exclusive=0),
            Form.FloatInterval(
                't_low', 't_high', 'divtime interval',
                '0', '5', low_inclusive=0, low_width_exclusive=0),
            Form.RadioGroup('legend_placement', 'plot legend location', [
                Form.RadioItem(TOPLEFT, 'top left'),
                Form.RadioItem(BOTTOMLEFT, 'bottom left'),
                Form.RadioItem(TOPRIGHT, 'top right', True),
                Form.RadioItem(BOTTOMRIGHT, 'bottom right')]),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('plot')

def compute_criterion(nstates, mu, t):
    """
    @param nstates: number of states in the jc69 process
    @param mu: the randomization rate of the process
    @param t: divergence time
    @return: the value of the criterion at the given divergence time
    """
    return math.log(1.0 + (nstates-1.0)*math.exp(-2*mu*t))

def make_table(args):
    """
    Make outputs to pass to RUtil.get_table_string.
    @param args: user args
    @return: matrix, headers
    """
    # define some variables
    t_low = args.t_low
    t_high = args.t_high
    if t_high <= t_low:
        raise ValueError('low time must be smaller than high time')
    ntimes = 100
    incr = (t_high - t_low) / (ntimes - 1)
    # define the numbers in the table
    arr = []
    for i in range(ntimes):
        t = t_low + i * incr
        a = compute_criterion(args.nstates_a, args.mu_a, t)
        b = compute_criterion(args.nstates_b, args.mu_b, t)
        row = [t, a, b]
        arr.append(row)
    headers = ['t', 'alpha', 'beta']
    return arr, headers

def get_response_content(fs):
    # legend labels
    label_a = 'N=%d mu=%f' % (fs.nstates_a, fs.mu_a)
    label_b = 'N=%d mu=%f' % (fs.nstates_b, fs.mu_b)
    arr, headers = make_table(fs)
    # compute the max value
    ymax = math.log(max(fs.nstates_a, fs.nstates_b))
    nfifths = int(math.floor(ymax * 5.0)) + 1
    ylim = RUtil.mk_call_str('c', 0, 0.2 * nfifths)
    # write the R script body
    out = StringIO()
    print >> out, RUtil.mk_call_str(
            'plot',
            'my.table$t',
            'my.table$alpha',
            type='"n"',
            ylim=ylim,
            xlab='"time"',
            ylab='"information"',
            main='"comparison of an information criterion for two processes"',
            )
    # draw some horizontal lines
    for i in range(nfifths+1):
        print >> out, RUtil.mk_call_str(
                'abline',
                h=0.2*i,
                col='"lightgray"',
                lty='"dotted"')
    colors = ('darkblue', 'darkred')
    for c, header in zip(colors, headers[1:]):
        print >> out, RUtil.mk_call_str(
                'lines',
                'my.table$t',
                'my.table$%s' % header,
                col='"%s"' % c,
                )
    legend_names = (label_a, label_b)
    legend_name_str = 'c(' + ', '.join('"%s"' % s for s in legend_names) + ')'
    legend_col_str = 'c(' + ', '.join('"%s"' % s for s in colors) + ')'
    legend_lty_str = 'c(' + ', '.join('1' for s in colors) + ')'
    print >> out, RUtil.mk_call_str(
            'legend',
            '"%s"' % fs.legend_placement,
            legend_name_str,
            col=legend_col_str,
            lty=legend_lty_str,
            )
    script_body = out.getvalue()
    # create the R plot image
    table_string = RUtil.get_table_string(arr, headers)
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script_body, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

