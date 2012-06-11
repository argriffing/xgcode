"""
Utility functions for interfacing with R.

Note that to use tikz with R,
you should put a tikzMetricsDictionary path into your .Rprofile dotfile.
For example
options( tikzMetricsDictionary='/home/username/.tikzMetricsDictionary' )
Otherwise it takes a few seconds to generate a temporary dictionary
every time you want to make a tikz.
"""

from StringIO import StringIO
from collections import defaultdict
import os
import re
import unittest
import subprocess

import Util

# pillaged from Carbone.py but more liberal
g_header_pattern = r'^[\.a-zA-Z][\.a-zA-Z0-9]*$'

#TODO allow RTable flexibility to not have row headers

#FIXME
# This is a horrible hack.
# It is a result of the intersection of two bad things.
# First, there is no system-wide configuration file for the python scripts.
# Second, the servbio hpc does not have R on the path
# which is passed to the python script from the mobyle caller.
# Note that I used to try to use
# /usr/local/apps/R/xeon/2.9.0/bin/R
# but this version connected to an invalid readline library.
g_rlocations = (
        'R',
        '/usr/local/apps/R/em64t/R-2.11.1/bin/R')


# This is an example R code from
# http://www.texample.net/tikz/examples/tikzdevice-demo/
#
g_stub = r"""
# Normal distribution curve
x <- seq(-4.5,4.5,length.out=100)
y <- dnorm(x)
#
# Integration points
xi <- seq(-2,2,length.out=30)
yi <- dnorm(xi)
#
# plot the curve
plot(x,y,type='l',col='blue',ylab='$p(x)$',xlab='$x$')
# plot the panels
lines(xi,yi,type='s')
lines(range(xi),c(0,0))
lines(xi,yi,type='h')
#
#Add some equations as labels
title(main="$p(x)=\\frac{1}{\\sqrt{2\\pi}}e^{-\\frac{x^2}{2}}$")
int <- integrate(dnorm,min(xi),max(xi),subdivisions=length(xi))
text(2.8, 0.3, paste("\\small$\\displaystyle\\int_{", min(xi),
"}^{", max(xi), "}p(x)dx\\approx", round(int[['value']],3),
'$', sep=''))
""".strip()

g_devices = {'pdf', 'postscript', 'png', 'tikz'}

class RError(Exception):
    """
    This is an error running R.
    """
    pass

class RExecError(Exception):
    """
    This is an error finding an R to run.
    """
    pass

class RTableError(Exception):
    """
    This is an R table parsing error.
    """
    pass

# pillaged from Carbone.py
def is_valid_header(h):
    return re.match(g_header_pattern, h)

# pillaged from Carbone.py
def validate_headers(headers):
    for h in headers:
        if not is_valid_header(h):
            raise RTableError('invalid column header: %s' % h)


def mk_call_str(name, *args, **kwargs):
    args_v = [str(v) for v in args]
    kwargs_v = ['%s=%s' % kv for kv in kwargs.items()]
    arr = args_v + kwargs_v
    return '%s(%s)' % (name, ', '.join(arr))

def run(pathname, rlocations=g_rlocations):
    """
    Run the R script.
    Redirect .Rout to stderr.
    The return code is 0 for success and 1 for failure.
    The returned stdout and stderr are strings not files.
    @param pathname: name of the R script
    @param rlocations: place to look for R
    @return: (returncode, r_stdout, r_stderr)
    """
    proc = None
    for rlocation in rlocations:
        # Note that we do not use --vanilla because we need ~/.Rprofile
        # to have the path to the tikz metrics dictionary
        # so that the tikz device does rebuild the dictionary each time.
        cmd = [
                rlocation, 'CMD', 'BATCH',
                '--no-save', '--no-restore',
                # the following flags would be added by --vanilla:
                #'--no-environ', '--no-site-file', '--no-init-file',
                '--silent', '--slave',
                pathname, '/dev/stderr']
        try:
            proc = subprocess.Popen(cmd,
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except OSError as e:
            continue
    if proc is None:
        raise RExecError('could not find R')
    proc_stdout, proc_stderr = proc.communicate()
    return proc.returncode, proc_stdout, proc_stderr

def run_with_table(table, user_data, callback):
    """
    @param table: the table string
    @param user_data: typically a fieldstorage-like object
    @param callback: this callback is like f(user_data, table_filename)
    @return: the R output as a string
    """
    retcode, r_out, r_err = run_with_table_verbose(
            table, user_data, callback)
    if retcode:
        raise RError(r_err)
    return r_err

def run_with_table_verbose(table, user_data, callback):
    """
    @param table: the table string
    @param user_data: typically a fieldstorage-like object
    @param callback: this callback is like f(user_data, table_filename)
    @return: returncode, r_stdout, r_stderr
    """
    # Create a temporary data table file for R.
    f_temp_table = Util.create_tmp_file(table)
    # Create a temporary R script file.
    script_content = callback(user_data, f_temp_table)
    f_temp_script = Util.create_tmp_file(script_content)
    # Call R.
    retcode, r_out, r_err = run(f_temp_script)
    # To facilitate debugging, only delete temporary files if R was successful.
    if not retcode:
        # Delete the temporary data table file.
        os.unlink(f_temp_table)
        # Delete the temporary script file.
        os.unlink(f_temp_script)
    # Return the R results.
    return retcode, r_out, r_err

def _get_device_specific_call(temp_plot_name, device_name,
        width=None, height=None):
    if device_name == 'tikz':
        call_string = 'system.time(tikz("%s"' % temp_plot_name
        if width:
            call_string += ', width=%f' % width
        if height:
            call_string += ', height=%f' % height
        call_string += '))'
    else:
        call_string = '%s("%s")' % (device_name, temp_plot_name)
    return call_string

def run_plotter_concise(table, user_script_content, device_name,
        width=None, height=None, keep_intermediate=False):
    """
    Raise an error instead of returning retcode and stderr and stdout.
    This should probably be the usual invocation and run_plotter
    should be changed to run_plotter_verbose, but this is a transitional form.
    """
    retcode, r_out, r_err, image_data = run_plotter(
            table, user_script_content, device_name,
            width, height, keep_intermediate)
    if retcode:
        raise RError(r_err)
    return image_data

def run_plotter(table, user_script_content, device_name,
        width=None, height=None, keep_intermediate=False):
    """
    The header and footer of the script are automatically included.
    The header reads the table and loads device-specific libraries
    (e.g. for tikz) and turns on the device.
    The footer turns off the device.
    @param table: the table string
    @param user_script_content: script without header or footer
    @param device_name: an R device function name
    @param width: optional width passed to tikz
    @param height: optional height passed to tikz
    @param keep_intermediate: a flag to keep the intermediate files
    @return: returncode, r_stdout, r_stderr, image_data
    """
    temp_table_name = Util.create_tmp_file(table, suffix='.table')
    temp_plot_name = Util.get_tmp_filename()
    s = StringIO()
    print >> s, 'my.table <- read.table("%s")' % temp_table_name
    if device_name == 'tikz':
        print >> s, 'require(tikzDevice)'
    print >> s, _get_device_specific_call(temp_plot_name, device_name,
            width, height)
    print >> s, user_script_content
    print >> s, 'dev.off()'
    script_content = s.getvalue()
    temp_script_name = Util.create_tmp_file(script_content, suffix='.R')
    retcode, r_out, r_err = run(temp_script_name)
    image_data = None
    if not retcode:
        if not keep_intermediate:
            os.unlink(temp_table_name)
            os.unlink(temp_script_name)
        try:
            with open(temp_plot_name, 'rb') as fin:
                image_data = fin.read()
        except IOError as e:
            raise RError(
                    'could not open the plot image file '
                    'that R was supposed to write')
        if not keep_intermediate:
            os.unlink(temp_plot_name)
    return retcode, r_out, r_err, image_data

def run_plotter_multiple_scripts(table, scripts, device_name,
        width=None, height=None):
    """
    @param table: the table string
    @param scripts: user scripts without header or footer
    @param device_name: an R device function name
    @param width: optional width passed to tikz
    @param height: optional height passed to tikz
    @return: returncode, r_stdout, r_stderr, image_data
    """
    #TODO: reorganize the code to combine this with run_plotter
    temp_table_name = Util.create_tmp_file(table)
    temp_plot_names = [Util.get_tmp_filename() for x in scripts]
    s = StringIO()
    print >> s, 'my.table <- read.table("%s")' % temp_table_name
    if device_name == 'tikz':
        print >> s, 'require(tikzDevice)'
    for (plot_name, script) in zip(temp_plot_names, scripts):
        print >> s, _get_device_specific_call(plot_name, device_name,
                width, height)
        print >> s, script
        print >> s, 'dev.off()'
    script_content = s.getvalue()
    temp_script_name = Util.create_tmp_file(script_content)
    retcode, r_out, r_err = run(temp_script_name)
    image_data_list = []
    if not retcode:
        os.unlink(temp_table_name)
        os.unlink(temp_script_name)
        try:
            for temp_plot_name in temp_plot_names:
                with open(temp_plot_name, 'rb') as fin:
                    image_data_list.append(fin.read())
        except IOError as e:
            raise RError(
                    'could not open the plot image file '
                    'that R was supposed to write')
        os.unlink(temp_plot_name)
    return retcode, r_out, r_err, image_data_list

def run_plotter_no_table(user_script_content, device_name,
        width=None, height=None):
    """
    @param user_script_content: script without header or footer
    @param device_name: an R device function name
    @param width: optional width passed to tikz
    @param height: optional height passed to tikz
    @return: returncode, r_stdout, r_stderr, image_data
    """
    temp_plot_name = Util.get_tmp_filename()
    s = StringIO()
    print >> s, _get_device_specific_preamble(temp_plot_name, device_name,
            width, height)
    print >> s, user_script_content
    print >> s, 'dev.off()'
    script_content = s.getvalue()
    temp_script_name = Util.create_tmp_file(script_content)
    retcode, r_out, r_err = run(temp_script_name)
    image_data = None
    if not retcode:
        os.unlink(temp_script_name)
        try:
            with open(temp_plot_name, 'rb') as fin:
                image_data = fin.read()
        except IOError as e:
            raise RError(
                    'could not open the plot image file '
                    'that R was supposed to write')
        os.unlink(temp_plot_name)
    return retcode, r_out, r_err, image_data

def get_table_string(M, column_headers):
    """
    Convert a row major rate matrix to a string representing an R table.
    @param M: a row major matrix
    @param column_headers: the labels of the data columns
    """
    if len(set(len(row) for row in M)) != 1:
        raise ValueError('all rows should have the same length')
    if len(M[0]) != len(column_headers):
        raise ValueError(
                'the number of columns does not match the number of headers')
    for header in column_headers:
        if '_' in header:
            raise ValueError(
                    'the header "%s" is invalid '
                    'because it has an underscore' % header)
    # define each line in the output string
    lines = []
    lines.append('\t'.join([''] + list(column_headers)))
    for i, row in enumerate(M):
        R_row = [float_to_R(value) for value in row]
        lines.append('\t'.join([str(i+1)] + R_row))
    return '\n'.join(lines)

def float_to_R(value):
    """
    Convert a python floating point value to a string usable by R.
    @param value: a floating point number
    @return: the R string representing the floating point number
    """
    if value == float('inf'):
        return 'Inf'
    else:
        return str(value)

def matrix_to_R_string(M):
    """
    @param M: a list of lists
    @return: a single line string suitable for copying and pasting into R
    """
    arr = []
    for row in M:
        for element in row:
            arr.append(str(element))
    return mk_call_str('matrix',
            mk_call_str('c', ', '.join(arr)),
            len(M), len(M[0]), byrow='TRUE')


class RTable:
    """
    Do things with a user-supplied R table in a way that expects errors.
    When the user inevitably wants to do something that causes an error,
    try to provide enough information so that they can decide
    how to change the table or their query.
    For example they may ask by name for a column that is not in the table,
    or multiple columns may have the same name,
    or they may want to use a column as a primary key
    but the column may have repeated values.
    This was initially implemented for some work with Ignazio Carbone
    who wanted access to R from the internet,
    specifically through Mobyle running on a Rocks cluster on CentOS.
    """
    def __init__(self, raw_lines):
        """
        Initialize some member variables.
        """
        # a list of header strings
        self.headers = []
        # a row major list of string data elements
        self.data = []
        # maps a header string to a column index
        self.h_to_i = {}
        # initialize the member variables.
        self._parse_r_table(raw_lines)
        self.h_to_i = dict((h, i+1) for i, h in enumerate(self.headers))

    def _parse_r_table(self, raw_lines):
        """
        Parse an R table into a header row and data rows.
        Each element of the data table is a string.
        """
        lines = Util.get_stripped_lines(raw_lines)
        header_line, data_lines = lines[0], lines[1:]
        self.headers = header_line.split()
        self.data = [line.split() for line in data_lines]
        nheaders = len(self.headers)
        if len(set(self.headers)) != nheaders:
            raise RTableError(
                    'multiple columns are labeled with the same header')
        for row in self.data:
            if len(row) != nheaders+1:
                raise RTableError(
                        'the header row has %d elements '
                        'and a data row has %d elements; '
                        'all data rows should have one more element '
                        'than the header row' % (nheaders, len(row)))
        validate_headers(self.headers)
    
    def header_to_column_index(self, header):
        if header not in self.h_to_i:
            raise RTableError(
                    'the column header %s was not found in the table' % header)
        return self.h_to_i[header]

    def header_to_column(self, header):
        column_index = self.header_to_column_index(header)
        column = [row[column_index] for row in self.data]
        return column

    def header_to_primary_column(self, header):
        """
        This is a strict variant of header_to_column.
        It requires the column entries to be unique.
        """
        column = self.header_to_column(header)
        d = defaultdict(int)
        for x in column:
            d[x] += 1
        repeated_keys = [k for k, v in d.items() if v > 1]
        if len(repeated_keys) > 5:
            raise RTableError(
                    '%d repeated keys '
                    'found in the primary column' % len(repeated_keys))
        elif repeated_keys:
            raise RTableError(
                'repeated keys '
                'in the primary column: %s' % ', '.join(repeated_keys))
        return column


class TestRUtil(unittest.TestCase):

    def test_mk_call_str(self):
        observed = mk_call_str('wat', 'x', 'y', z='foo', w='bar')
        expected = 'wat(x, y, z=foo, w=bar)'
        self.assertEqual(observed, expected)

    def test_matrix_to_R_string(self):
        observed = matrix_to_R_string([[1, 2, 3], [4, 5, 6]])
        expected = 'matrix(c(1, 2, 3, 4, 5, 6), 2, 3, byrow=TRUE)'
        self.assertEquals(expected, observed)

    def test_exec_error(self):
        self.assertRaises(RExecError, run, 'whatever.R', 'bad_r_location')


if __name__ == '__main__':
    unittest.main()

