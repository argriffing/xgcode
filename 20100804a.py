"""Scatter plot 2D given an R table with two categorical variables.
"""

from StringIO import StringIO
import os
import tempfile
import colorsys

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import RUtil
import Carbone
import iterutils
import const

g_tags = ['pca:plot']

g_default = const.read('20100709a')

g_colorbrewer_set1 = [
    "#E41A1C", "#377EB8", "#4DAF4A",
    "#984EA3", "#FF7F00", "#FFFF33",
    "#A65628", "#F781BF", "#999999"]

#FIXME multiple output types

def rgb_floats_to_ints(rgb_floats):
    return tuple(int(round(x*255)) for x in rgb_floats)

def create_R_palette(n):
    """
    @param n: make a list of this many colors
    @return: a list of R colors like #E41A1C
    """
    if n < 10:
        return g_colorbrewer_set1[:n]
    increment = 1.0 / (1 + n)
    hues = [i*increment for i in range(n)]
    # get rgb values as triples of floats between 0 and 1
    rgbs = [colorsys.hsv_to_rgb(h, 1.0, 1.0) for h in hues]
    # get rgb values as triples of ints between 0 and 255
    rgbs = [rgb_floats_to_ints(rgb) for rgb in rgbs]
    # get rgb values as strings
    rgbs = ['#%02X%02X%02X' % rgb for rgb in rgbs]
    return rgbs

class MissingError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default),
            Form.Sequence('axes',
                'numerical variables defining axes of the scatter plot',
                ('pc1', 'pc2')),
            Form.SingleLine('shape',
                'categorical variable defining the shape of a dot',
                'species'),
            Form.SingleLine('color',
                'categorical variable defining the color of a dot',
                'location'),
            Form.Float('size', 'size of a plotted point', '1.5'),
            Form.Sequence('symbol_legend_pos',
                'position of the symbol legend',
                ('-10', '-1')),
            Form.Sequence('color_legend_pos',
                'position of the color legend',
                ('0', '-1')),
            Form.ImageFormat(),
            #Form.RadioGroup('out_type', 'output type', [
                #Form.RadioItem('show_image', 'image', True),
                #Form.RadioItem('show_table', 'R table'),
                #Form.RadioItem('show_script', 'R script')]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('%s.%s.pca.2d', ['shape', 'color'])

def get_response_content(fs):
    # create a response that depends on the requested output type
    #if fs.show_image:
    if 1:
        content = process(fs, fs.table.splitlines())
        #ext = Form.g_imageformat_to_ext[fs.imageformat]
        #filename = '.'.join((fs.shape, fs.color, 'pca', '2d', ext))
        #contenttype = Form.g_imageformat_to_contenttype[fs.imageformat]
    """
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
    """
    return content

class ColorInfo:
    def __init__(self, header, column):
        """
        @param header: the name of the variable to be represented by color
        @param column: a list of categorical values
        """
        self.header = header
        self.values = column[:]
        self.unique_values = list(iterutils.unique_everseen(column))
        self.unique_colors = create_R_palette(len(self.unique_values))
        value_to_color = dict(zip(self.unique_values, self.unique_colors))
        self.colors = [value_to_color[v] for v in column]
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
        return 'c(' + ','.join('"%s"' % x for x in self.unique_colors) + ')'

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
        self.axis_headers = args.axes
        # verify the number of axis headers
        if len(self.axis_headers) != 2:
            raise ValueError('expected two axis column headers')
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
                axis_list = Carbone.get_numeric_column(data, index)
            except Carbone.NumericError:
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
        nrows = len(self.color_info.colors)
        header_row = ['x', 'y', 'color', 'symbol']
        data_rows = zip(
                range(1, nrows+1),
                self.axis_lists[0],
                self.axis_lists[1],
                ['"%s"' % x for x in self.color_info.colors],
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
            symbol_legend_pos = Util.get_coordinate_pair(args.symbol_legend_pos)
        except Util.CoordinatePairError as e:
            msg = 'symbol legend position error: ' + str(e)
            raise ValueError(msg)
        # get the color legend location
        try:
            color_legend_pos = Util.get_coordinate_pair(args.color_legend_pos)
        except Util.CoordinatePairError as e:
            msg = 'color legend position error: ' + str(e)
            raise ValueError(msg)
        # get the image function
        image_function = Form.g_imageformat_to_r_function[args.imageformat]
        # get the R codes
        rcodes = [
            "mytable <- read.table('%s')" % temp_table_filename,
            "%s('%s')" % (image_function, temp_plot_filename),
            # create the scatterplot
            "myplot <- plot(mytable$x, mytable$y, mar=c(5,3,4,3),",
            "xlab = '%s'," % self.axis_headers[0],
            "ylab = '%s'," % self.axis_headers[1],
            "type = 'p', pch = ' ')",
            # define symbols colors and sizes
            "points(mytable$x, mytable$y,",
            "pch=mytable$symbol,"
            'bg=as.vector(mytable$color),',
            'col=as.vector(mytable$color),',
            "cex=%s)" % args.size,
            # symbol legend
            "legend(%s, %s," % symbol_legend_pos,
            'pch=' + self.shape_info.get_legend_pch() + ',',
            "yjust = 0,",
            'legend=' + self.shape_info.get_legend_legend() + ',',
            'title="%s",' % self.shape_info.header,
            "cex = 1.1)",
            # color legend
            "legend(%s, %s," % color_legend_pos,
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
    temp_plot_name = Util.get_tmp_filename()
    # Create a temporary R script file.
    f_temp_script = tempfile.NamedTemporaryFile(delete=False)
    script = plot_info.get_script(args, temp_plot_name, f_temp_table.name)
    print >> f_temp_script, script
    f_temp_script.close()
    # Call R.
    retcode, r_out, r_err = RUtil.run(f_temp_script.name)
    if retcode:
        raise ValueError('R error:\n' + r_err)
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
