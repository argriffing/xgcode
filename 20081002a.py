"""Sample some points from interlocking spirals in two dimensions.

Each output row gives the x, y, and possibly the label associated with a point.
The label is either 1 or -1 depending on the group to which the point belongs.
This is not so useful.
"""

from StringIO import StringIO
import math

from SnippetUtil import HandlingError
import RUtil
import SpiralSampler
import Form
import FormOut

#FIXME clarify output format

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = [
            Form.Integer('npoints', 'sample this many points per group',
                100, low=0, high=10000),
            Form.Float('stddev', 'standard deviation of the noise',
                0.1, low_exclusive=0),
            Form.CheckGroup('label_options', 'label options', [
                Form.CheckItem('add_labels',
                    'add group label to each row', True)]),
            Form.RadioGroup('format', 'format options', [
                Form.RadioItem('raw', 'rows of tab separated values', True),
                Form.RadioItem('table', 'R table format')]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # unpack some options
    npoints = fs.npoints
    stddev = fs.stddev
    # define the data rows and the headers
    if fs.add_labels:
        headers = ('x', 'y', 'label')
        data_rows = list(SpiralSampler.gen_labeled_points(npoints, stddev))
    else:
        headers = ('x', 'y')
        data_rows = list(SpiralSampler.gen_points(npoints, stddev))
    # begin the response
    if fs.raw:
        lines = []
        for data_row in data_rows:
            line = '\t'.join(str(x) for x in data_row)
            lines.append(line)
        response_text = '\n'.join(lines)
    elif fs.table:
        response_text = RUtil.get_table_string(data_rows, headers)
    # return the response
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'spiral.table')
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
    return response_headers, response_text

