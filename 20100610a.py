"""Scatter plot 3D given an R table with one categorical and one numerical var.
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
                'numerical variable defining the color of a dot',
                'temperature'),
            Form.Float('size', 'size of a plotted point', '1.5'),
            Form.SingleLine('legend_pos', 'position of the symbol legend',
                '0 0 0'),
            Form.CheckGroup('pixel_options', 'more plotting details', [
                Form.CheckItem('endpoint_ticks',
                    'add ticks to the endpoints of the colorbar', True)]),
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


class PlotInfo:
    def __init__(self, args, headers, data):
        """
        @param args: user args from web or cmdline
        @param data: 
        """
        # map the column header to the column index
        self.h_to_i = dict((h, i+1) for i, h in enumerate(headers))
        # init the info
        self._init_axes(args, headers, data)
        self._init_colors(args, headers, data)
        self._init_shapes(args, headers, data)
        self._init_unique_shapes()

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

    def _init_colors(self, args, headers, data):
        """
        Colors are numeric, and use whatever gradient is built into R.
        """
        self.color_header = args.color
        if self.color_header not in headers:
            msg = 'bad color column header: ' + self.color_header
            raise ValueError(msg)
        index = self.h_to_i[self.color_header]
        try:
            self.color_list = get_numeric_column(data, index)
        except NumericError:
            msg_a = 'expected the color column %s ' % self.color_header
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)

    def _init_shapes(self, args, headers, data):
        """
        Shapes are categorical.
        """
        self.shape_header = args.shape
        if self.shape_header not in headers:
            msg = 'bad shape column header: ' + self.shape_header
            raise ValueError(msg)
        index = self.h_to_i[self.shape_header]
        self.shape_list = zip(*data)[index]

    def _init_unique_shapes(self):
        self.unique_shapes = list(iterutils.unique_everseen(self.shape_list))

    def get_augmented_table_lines(self):
        """
        This is given to R.
        """
        nrows = len(self.shape_list)
        header_row = ['x', 'y', 'z', 'color', 'symbol']
        data_rows = zip(
                range(1, nrows+1),
                self.axis_lists[0],
                self.axis_lists[1],
                self.axis_lists[2],
                self.color_list,
                [self.unique_shapes.index(x) for x in self.shape_list])
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
        legend_pos = tuple(float(x) for x in args.legend_pos.split())
        if len(legend_pos) != 3:
            raise ValueError('invalid legend position: %s' % args.legend_pos)
        # get the unique locations and species
        symbol_legend_string = ', '.join("'%s'" % x for x in self.unique_shapes)
        color_legend_string = self.color_header
        # add color legend endpoint axis
        if args.endpoint_ticks:
            s = 'mytable$color'
            color_axis = 'axis(1, c(axTicks(1), min(%s), max(%s)))' % (s, s)
        else:
            color_axis = 'axis(1)'
        rcodes = [
            "require('scatterplot3d')",
            "mytable <- read.table('%s')" % temp_table_filename,
            "%s('%s')" % (args.imageformat, temp_plot_filename),
            # rename some variables for compatibility with the template
            "Year <- mytable$x",
            "Latitude <- mytable$y",
            "Risk <- mytable$z",
            "Prec <- mytable$color",
            # stack two plots vertically
            "layout(cbind(1:2, 1:2), heights = c(7, 1))",
            # create the color gradient
            "prc <- hsv((prc <- 0.7*Prec/diff(range(Prec))) - min(prc) + 0.3)",
            # create the scatterplot
            "s3d <- scatterplot3d(Year, Latitude, Risk, mar = c(5, 3, 4, 3),",
            "xlab = '%s'," % self.axis_headers[0],
            "ylab = '%s'," % self.axis_headers[1],
            "zlab = '%s'," % self.axis_headers[2],
            "type = 'p', pch = ' ')",
            # define symbols colors and sizes
            "s3d$points(Year, Latitude, Risk,",
            "pch=mytable$symbol, bg=prc, col=prc, cex=%s)" % args.size,
            # define x y and z as Year, Latitude and Risk
            "s3d.coords <- s3d$xyz.convert(Year, Latitude, Risk)",
            # symbol legend
            "legend(s3d$xyz.convert(%s, %s, %s)," % legend_pos,
            "pch=0:%s, yjust = 0," % (len(symbol_legend_string)-1),
            "legend=c(%s)," % symbol_legend_string,
            "cex = 1.1)",
            # set margins
            "par(mar=c(5, 3, 0, 3))",
            # draw the plot
            "plot(seq(min(Prec), max(Prec), length = 100),",
            "rep(0, 100), pch = 15,",
            "axes = FALSE,",
            "xlab = '%s'," % color_legend_string,
            "ylab = '', col = hsv(seq(0.3, 1, length = 100)))",
            # draw the axis onto the color legend
            color_axis,
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

def main(args):
    with open(os.path.expanduser(args.table)) as fin:
        sys.stdout.write(process(args, fin))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', required=True,
            help='an R table with a specific format')
    parser.add_argument('--imageformat', default='png',
            choices=['png', 'pdf', 'postscript', 'svg'],
            help='output image file format')
    parser.add_argument('--shape', required=True,
            help='categorical variable defining the shape of a dot')
    parser.add_argument('--color', required=True,
            help='categorical variable defining the color of a dot')
    parser.add_argument('--size', default=1.5, type=float,
            help='the size of each plotted point')
    parser.add_argument('--legend_pos', default='0 0 0',
            help='position of the symbol legend')
    parser.add_argument('--axes', default='pc1 pc2 pc3',
            help='numerical variables defining axes of the 3D plot')
    main(parser.parse_args())

