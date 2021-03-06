"""
Scatter plot 3D given an R table with one categorical and one numerical var.
"""

from StringIO import StringIO
import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import RUtil
import iterutils
import const

g_tags = ['pca:plot']

g_default = const.read('20100709a')

#FIXME multiple output types

class MissingError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default),
            Form.Sequence('axes',
                'numerical variables defining axes of the 3D plot',
                ('pc1', 'pc2', 'pc3')),
            Form.SingleLine('shape',
                'categorical variable defining the shape of a dot',
                'species'),
            Form.SingleLine('color',
                'numerical variable defining the color of a dot',
                'temperature'),
            Form.Float('size', 'size of a plotted point', '1.5'),
            Form.Sequence('legend_pos', 'position of the symbol legend',
                ('0', '0', '0')),
            Form.CheckGroup('pixel_options', 'more plotting details', [
                Form.CheckItem('endpoint_ticks',
                    'add ticks to the endpoints of the colorbar', True)]),
            Form.ImageFormat()]
            #Form.RadioGroup('out_type', 'output type', [
                #Form.RadioItem('show_image', 'image', True),
                #Form.RadioItem('show_table', 'R table'),
                #Form.RadioItem('show_script', 'R script')]),
    return form_objects

def get_form_out():
    return FormOut.Image('%s.%s.pca.3d', ['shape', 'color'])

def get_response_content(fs):
    # create a response that depends on the requested output type
    #if fs.show_image:
    if 1:
        content = process(fs, fs.table.splitlines())
        ext = Form.g_imageformat_to_ext[fs.imageformat]
        filename = '.'.join((fs.shape, fs.color, 'pca', '3d', ext))
        contenttype = Form.g_imageformat_to_contenttype[fs.imageformat]
    """
    else:
        # read the table
        rtable = RUtil.RTable(fs.table.splitlines())
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
        self.axis_headers = args.axes
        # verify the number of axis headers
        if len(self.axis_headers) != 3:
            raise ValueError('expected three axis column headers')
        # verify the axis header contents
        bad_axis_headers = set(self.axis_headers) - set(headers)
        if bad_axis_headers:
            raise ValueError(
                    'bad axis column headers: ' + ', '.join(bad_axis_headers))
        self.axis_lists = []
        for h in self.axis_headers:
            index = self.h_to_i[h]
            try:
                axis_list = Carbone.get_numeric_column(data, index)
            except Carbone.NumericError:
                raise ValueError(
                        'expected the axis column %s '
                        'to be numeric' % h)
            self.axis_lists.append(axis_list)

    def _init_colors(self, args, headers, data):
        """
        Colors are numeric, and use whatever gradient is built into R.
        """
        self.color_header = args.color
        if self.color_header not in headers:
            raise ValueError('bad color column header: ' + self.color_header)
        index = self.h_to_i[self.color_header]
        try:
            self.color_list = Carbone.get_numeric_column(data, index)
        except Carbone.NumericError:
            raise ValueError(
                    'expected the color column %s '
                    'to be numeric' % self.color_header)

    def _init_shapes(self, args, headers, data):
        """
        Shapes are categorical.
        """
        self.shape_header = args.shape
        if self.shape_header not in headers:
            raise ValueError('bad shape column header: ' + self.shape_header)
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


    def get_script(self, args):
        """
        @param args: from cmdline or web
        """
        # get the symbol legend location
        try:
            legend_pos = Util.get_coordinate_triple(args.legend_pos)
        except Util.CoordinateTripleError as e:
            raise ValueError('legend position error: ' + str(e))
        # get the unique locations and species
        symbol_legend_string = ', '.join("'%s'" % x for x in self.unique_shapes)
        color_legend_string = self.color_header
        # add color legend endpoint axis
        if args.endpoint_ticks:
            s = 'my.table$color'
            color_axis = 'axis(1, c(axTicks(1), min(%s), max(%s)))' % (s, s)
        else:
            color_axis = 'axis(1)'
        # get the image function
        rcodes = [
            "require('scatterplot3d')",
            # rename some variables for compatibility with the template
            "Year <- my.table$x",
            "Latitude <- my.table$y",
            "Risk <- my.table$z",
            "Prec <- my.table$color",
            # stack two plots vertically
            "layout(cbind(1:2, 1:2), heights = c(7, 1.5))",
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
            "pch=my.table$symbol, bg=prc, col=prc, cex=%s)" % args.size,
            # define x y and z as Year, Latitude and Risk
            "s3d.coords <- s3d$xyz.convert(Year, Latitude, Risk)",
            # symbol legend
            "legend(s3d$xyz.convert(%s, %s, %s)," % legend_pos,
            "pch=0:%s, yjust = 0," % (len(symbol_legend_string)-1),
            "legend=c(%s)," % symbol_legend_string,
            "cex = 1.1)",
            # set margins
            "par(mar=c(5, 3, 0.5, 3))",
            # draw the plot
            "plot(seq(min(Prec), max(Prec), length = 100),",
            "rep(0, 100), pch = 15, cex = 2,",
            "axes = FALSE,",
            "xlab = '%s'," % color_legend_string,
            "ylab = '', col = hsv(seq(0.3, 1, length = 100)))",
            # draw the axis onto the color legend
            color_axis,
            ]
        return '\n'.join(rcodes)


def process(args, table_lines):
    """
    @param args: command line or web input
    @param table_lines: input lines
    @return: the image data as a string
    """
    # get the table string
    rtable = RUtil.RTable(table_lines)
    plot_info = PlotInfo(args, rtable.headers, rtable.data)
    augmented_lines = plot_info.get_augmented_table_lines()
    table_string = '\n'.join(augmented_lines)
    # get the script string
    script_string = plot_info.get_script(args)
    # get the device
    device = Form.g_imageformat_to_r_function[args.imageformat]
    # run R and get the image data
    image_data = RUtil.run_plotter_concise(table_string, script_string, device)
    # return the image data
    return image_data

