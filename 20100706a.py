"""Scatter plot 3D given an R table with two categorical variables.
"""

from StringIO import StringIO
import os
import subprocess
from subprocess import PIPE
import tempfile

import argparse

from SnippetUtil import HandlingError
import Form
import Util
import Carbone
import iterutils

g_default_rows = [
        ['otu', 'species', 'location', 'temperature', 'precipitation',
            'pc1', 'pc2', 'pc3'],
        [1, 'IC100', 'Ap', 'GA', 15.0, 600.0,
            -2.8053476259, 0.556532380058, -6.17891756957],
        [2, 'IC101', 'Ap', 'GA', 15.0, 600.0,
            -2.8053476259, 0.556532380058, -6.17891756956],
        [3, 'IC102', 'Ap', 'GA', 15.0, 600.0,
            -2.80534762591, 0.556532380059, -6.17891756957],
        [455, 'IC577', 'Ac', 'AR', 25.0, 400.0,
            -13.7544736082, -7.16259232881, 7.0902951321],
        [456, 'IC580', 'Ac', 'AR', 25.0, 400.0,
            3.56768959361, 0.385873934264, 1.23981735331],
        [457, 'IC591', 'Ac', 'AR', 25.0, 400.0,
            -11.6455270418, -5.710582374, 5.60835091179]]

g_default_lines = ['\t'.join(str(x) for x in row) for row in g_default_rows]
g_default_string = '\n'.join(g_default_lines)

def my_mktemp():
    """
    The old tempfile.mktemp is deprecated.
    This function will create the name of a file
    that an external process can create.
    @return: the name of a nonexistent temporary file
    """
    f = tempfile.NamedTemporaryFile(delete=False)
    name = f.name
    f.close()
    os.unlink(name)
    return name

class MissingError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default_string),
            Form.SingleLine('axes',
                'numerical variables defining axes of the 3D plot',
                ' '.join(('pc1', 'pc2', 'pc3'))),
            Form.SingleLine('shape',
                'categorical variable defining the shape of a dot',
                'species'),
            Form.SingleLine('color',
                'categorical variable defining the color of a dot',
                'location'),
            Form.Float('size', 'size of a plotted point', '1.5'),
            Form.SingleLine('symbol_legend_pos',
                'position of the symbol legend',
                '0 0 0'),
            Form.SingleLine('color_legend_pos',
                'position of the color legend',
                '0 0 0'),
            Form.ImageFormat(),
            Form.RadioGroup('out_type', 'output type', [
                Form.RadioItem('show_image', 'image', True),
                Form.RadioItem('show_table', 'R table'),
                Form.RadioItem('show_script', 'R script')]),
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # create a response that depends on the requested output type
    if fs.show_image:
        content = process(fs, fs.table.splitlines())
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        filename = '.'.join((fs.shape, fs.color, 'pca', '3d', ext))
        contenttype = Form.g_imageformat_to_contenttype[fs.imageformat]
    else:
        # read the table
        rtable = Carbone.RTable(fs.table.splitlines())
        header_row = rtable.headers
        data_rows = rtable.data
        # Do a more stringent check of the column headers.
        for h in header_row:
            if not Carbone.is_valid_header(h):
                msg = 'invalid column header: %s' % h
                raise ValueError(msg)
        plot_info = PlotInfo(fs, header_row, data_rows)
        if fs.show_table:
            content = '\n'.join(plot_info.get_augmented_table_lines()) + '\n'
            contenttype = 'text/plain'
            filename = 'out.table'
        elif fs.show_script:
            stub_image_name = 'stub-image-filename.' + fs.imageformat
            stub_table_name = 'stub-table-filename.table'
            content = plot_info.get_script(
                    fs, stub_image_name, stub_table_name) + '\n'
            contenttype = 'text/plain'
            filename = 'script.R'
    # return the response
    disposition = '%s; filename=%s' % (fs.contentdisposition, filename)
    response_headers = [
            ('Content-Type', contenttype),
            ('Content-Disposition', disposition)]
    return response_headers, content

class NumericError(Exception): pass
class LegendPositionError(Exception): pass

def get_numeric_column(data, index):
    """
    @param data: row major list of lists of numbers as strings
    @param index: column index
    @return: list of floats
    """
    strings = zip(*data)[index]
    try:
        floats = [float(x) for x in strings]
    except ValueError, v:
        raise NumericError
    return floats

def get_hues(n):
    increment = 1.0 / (1 + n)
    return [i*increment for i in range(n)]

def get_legend_position(legend_position_string):
    """
    @param legend_position_string: something like '10 20 -3.2'
    @return: something like (10.0, 20.0, -3.2)
    """
    if not legend_position_string:
        msg = 'expected some coordinates'
        raise LegendPositionError(msg)
    triple = []
    for s in legend_position_string.split():
        try:
            x = float(s)
        except ValueError:
            msg = 'expected a number but got "%s"' % s
            raise LegendPositionError(msg)
        triple.append(x)
    if len(triple) != 3:
        msg = 'expected three whitespace separated coordinates'
        raise LegendPositionError(msg)
    return tuple(triple)

class ColorInfo:
    def __init__(self, header, column):
        """
        @param header: the name of the variable to be represented by color
        @param column: a list of categorical values
        """
        self.header = header
        self.values = column[:]
        self.unique_values = list(iterutils.unique_everseen(column))
        self.unique_hues = get_hues(len(self.unique_values))
        value_to_hue = dict(zip(self.unique_values, self.unique_hues))
        self.hues = [value_to_hue[v] for v in column]
        # pch fifteen is a solid block
        self.pch = 15
    def get_legend_pch(self):
        """
        Return a string to write into an R script.
        """
        return 'rep(%d, %d)' % (self.pch, len(self.unique_values))
    def get_legend_legend(self):
        """
        Return a string to write into an R script.
        @return: the legend parameter of the legend function for R
        """
        return 'c(' + ','.join('"%s"' % x for x in self.unique_values) + ')'
    def get_legend_col(self):
        """
        Return a string to write into an R script.
        """
        return 'hsv(c(%s))' % ','.join(str(x) for x in self.unique_hues)
    def get_dot_col(self):
        """
        Return a string to write into an R script.
        """
        return 'hsv(c(%s))' % ','.join(str(x) for x in self.hues)

class ShapeInfo:
    def __init__(self, header, column):
        """
        @param header: the name of the variable to be represented by shape
        @param column: a list of categorical values
        """
        self.header = header
        self.values = column[:]
        self.unique_values = list(iterutils.unique_everseen(column))
        # for now just use sequential integers starting at zero
        self.unique_pchs = range(len(self.unique_values))
        value_to_pch = dict(zip(self.unique_values, self.unique_pchs))
        self.pchs = [value_to_pch[v] for v in column]
    def get_legend_pch(self):
        """
        Return a string to write into an R script.
        @return: the pch parameter of the legend function for R
        """
        return 'c(' + ','.join(str(x) for x in self.unique_pchs) + ')'
    def get_legend_legend(self):
        """
        Return a string to write into an R script.
        @return: the legend parameter of the legend function for R
        """
        return 'c(' + ','.join('"%s"' % x for x in self.unique_values) + ')'

class PlotInfo:
    def __init__(self, args, headers, data):
        """
        @param args: user args from web or cmdline
        @param data: 
        """
        # map the column header to the column index
        self.h_to_i = dict((h, i+1) for i, h in enumerate(headers))
        # init the info about how to make the plot
        self._init_axes(args, headers, data)
        color_column = self._get_color_column(args.color, headers, data)
        shape_column = self._get_shape_column(args.shape, headers, data)
        self.color_info = ColorInfo(args.color, color_column)
        self.shape_info = ShapeInfo(args.shape, shape_column)

    def _init_axes(self, args, headers, data):
        # read the axes
        self.axis_headers = args.axes.split()
        # verify the number of axis headers
        if len(self.axis_headers) != 3:
            raise ValueError('expected three axis column headers')
        # verify the axis header contents
        bad_axis_headers = set(self.axis_headers) - set(headers)
        if bad_axis_headers:
            msg_a = 'bad axis column headers: '
            msg_b = ', '.join(bad_axis_headers)
            raise ValueError(msg_a + msg_b)
        self.axis_lists = []
        for h in self.axis_headers:
            index = self.h_to_i[h]
            try:
                axis_list = get_numeric_column(data, index)
            except NumericError:
                msg_a = 'expected the axis column %s ' % h
                msg_b = 'to be numeric'
                raise ValueError(msg_a + msg_b)
            self.axis_lists.append(axis_list)

    def _get_color_column(self, header, headers, data):
        """
        Colors are categorical.
        """
        if header not in headers:
            msg = 'bad color column header: ' + header
            raise ValueError(msg)
        index = self.h_to_i[header]
        return zip(*data)[index]

    def _get_shape_column(self, header, headers, data):
        """
        Shapes are categorical.
        """
        if header not in headers:
            msg = 'bad shape column header: ' + header
            raise ValueError(msg)
        index = self.h_to_i[header]
        return zip(*data)[index]

    def get_augmented_table_lines(self):
        """
        This is given to R.
        """
        nrows = len(self.color_info.hues)
        header_row = ['x', 'y', 'z', 'color', 'symbol']
        data_rows = zip(
                range(1, nrows+1),
                self.axis_lists[0],
                self.axis_lists[1],
                self.axis_lists[2],
                self.color_info.hues,
                self.shape_info.pchs)
        header_line = '\t'.join(str(x) for x in header_row)
        data_lines = ['\t'.join(str(x) for x in row) for row in data_rows]
        return [header_line] + data_lines

    def get_script(self, args, temp_plot_filename, temp_table_filename):
        """
        @param args: from cmdline or web
        @param temp_plot_name: a pathname
        @param temp_table_name: a pathname
        """
        # get the symbol legend location
        try:
            symbol_legend_pos = get_legend_position(args.symbol_legend_pos)
        except LegendPositionError as e:
            msg_a = 'symbol legend position error for '
            msg_b = '"%s" : %s' % (args.symbol_legend_pos, e)
            raise ValueError(msg_a + msg_b)
        # get the color legend location
        try:
            color_legend_pos = get_legend_position(args.color_legend_pos)
        except LegendPositionError as e:
            msg_a = 'color legend position error for '
            msg_b = '"%s" : %s' % (args.symbol_legend_pos, e)
            raise ValueError(msg_a + msg_b)
        # get the R codes
        rcodes = [
            "require('scatterplot3d')",
            "mytable <- read.table('%s')" % temp_table_filename,
            "%s('%s')" % (args.imageformat, temp_plot_filename),
            # rename some variables for compatibility with the template
            "Year <- mytable$x",
            "Latitude <- mytable$y",
            "Risk <- mytable$z",
            "Prec <- mytable$color",
            # create the scatterplot
            "s3d <- scatterplot3d(Year, Latitude, Risk, mar = c(5, 3, 4, 3),",
            "xlab = '%s'," % self.axis_headers[0],
            "ylab = '%s'," % self.axis_headers[1],
            "zlab = '%s'," % self.axis_headers[2],
            "type = 'p', pch = ' ')",
            # define symbols colors and sizes
            "s3d$points(Year, Latitude, Risk,",
            "pch=mytable$symbol,"
            'bg=' + self.color_info.get_dot_col() + ',',
            'col=' + self.color_info.get_dot_col() + ',',
            "cex=%s)" % args.size,
            # define x y and z as Year, Latitude and Risk
            "s3d.coords <- s3d$xyz.convert(Year, Latitude, Risk)",
            # symbol legend
            "legend(s3d$xyz.convert(%s, %s, %s)," % symbol_legend_pos,
            'pch=' + self.shape_info.get_legend_pch() + ',',
            "yjust = 0,",
            'legend=' + self.shape_info.get_legend_legend() + ',',
            'title="%s",' % self.shape_info.header,
            "cex = 1.1)",
            # color legend
            "legend(s3d$xyz.convert(%s, %s, %s)," % color_legend_pos,
            'pch=' + self.color_info.get_legend_pch() + ',',
            "yjust = 0,",
            'legend=' + self.color_info.get_legend_legend() + ',',
            'col=' + self.color_info.get_legend_col() + ',',
            'title="%s",' % self.color_info.header,
            "cex = 1.1)",
            # write the plot
            "dev.off()"]
        return '\n'.join(rcodes)


def process(args, table_lines):
    """
    @param args: command line or web input
    @param table_lines: input lines
    @return: the image data as a string
    """
    rtable = Carbone.RTable(table_lines)
    header_row = rtable.headers
    data_rows = rtable.data
    # Do a more stringent check of the column headers.
    for h in header_row:
        if not Carbone.is_valid_header(h):
            msg = 'invalid column header: %s' % h
            raise ValueError(msg)
    # Read the relevant columns and their labels.
    plot_info = PlotInfo(args, header_row, data_rows)
    # Get info for the temporary data
    augmented_lines = plot_info.get_augmented_table_lines()
    # Create a temporary data table file for R.
    f_temp_table = tempfile.NamedTemporaryFile(delete=False)
    print >> f_temp_table, '\n'.join(augmented_lines)
    f_temp_table.close()
    # Create a temporary pathname for the plot created by R.
    temp_plot_name = my_mktemp()
    # Create a temporary R script file.
    f_temp_script = tempfile.NamedTemporaryFile(delete=False)
    script = plot_info.get_script(args, temp_plot_name, f_temp_table.name)
    print >> f_temp_script, script
    f_temp_script.close()
    # Call R.
    cmd = [
            'R', 'CMD', 'BATCH',
            # turn off as much output as possible
            '--vanilla', '--slave', '--silent',
            f_temp_script.name,
            # avoid writing the .Rout file
            '/dev/null']
    proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE)
    r_out, r_err = proc.communicate()
    # Delete the temporary data table file.
    os.unlink(f_temp_table.name)
    # Delete the temporary script file.
    os.unlink(f_temp_script.name)
    # Read the image file.
    try:
        with open(temp_plot_name, 'rb') as fin:
            image_data = fin.read()
    except IOError, e:
        raise HandlingError('the R call seems to not have created the plot')
    # Delete the temporary image file.
    os.unlink(temp_plot_name)
    # Return the image data as a string.
    return image_data