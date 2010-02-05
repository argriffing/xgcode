"""Sample some points from interlocking spirals in two dimensions.

Each output row gives the x, y, and possibly the label associated with a point.
The label is either 1 or -1 depending on the group to which the point belongs.
This is not so useful.
"""

from StringIO import StringIO
import math

from SnippetUtil import HandlingError
import Form
import RUtil
import SpiralSampler

def get_form():
    """
    @return: a list of form objects
    """
    form_objects = [
            Form.Integer('npoints', 'sample this many points per group', 100, low=0, high=10000),
            Form.Float('stddev', 'standard deviation of the noise', 0.1, low_exclusive=0),
            Form.CheckGroup('label_options', 'label options', [
                Form.CheckItem('add_labels', 'add the group label to each row', True)]),
            Form.RadioGroup('format', 'format options', [
                Form.RadioItem('raw', 'rows of tab separated values', True),
                Form.RadioItem('table', 'R table format')]),
            Form.RadioGroup('contentdisposition', 'delivery options', [
                Form.RadioItem('inline', 'view', True),
                Form.RadioItem('attachment', 'download')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # define the data rows and the headers
    if fs.add_labels:
        headers = ('x', 'y', 'label')
        data_rows = list(SpiralSampler.gen_labeled_points(fs.npoints, fs.stddev))
    else:
        headers = ('x', 'y')
        data_rows = list(SpiralSampler.gen_points(fs.npoints, fs.stddev))
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
    response_headers = [('Content-Type', 'text/plain')]
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.contentdisposition, 'spiral.table')))
    return response_headers, response_text

