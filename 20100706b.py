"""Cluster using k-means. [UNFINISHED]
"""

from StringIO import StringIO
import os

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

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default_string),
            Form.SingleLine('axes', 'column labels of Euclidean axes',
                ' '.join(('pc1', 'pc2', 'pc3'))),
            Form.Integer('k', 'number of clusters',
                2, low=2),
            Form.SingleLine('annotation', 'header of added column',
                'cluster'),
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
        # init the info about how to make the plot
        self._init_axes(args, headers, data)

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
    # Return the image data as a string.
    return image_data
