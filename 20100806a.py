"""Given an R table, remove rows that are near-neighbors.
"""

from StringIO import StringIO

import numpy as np

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import const

g_tags = ['pca:misc']

g_default = const.read('20100709a')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default),
            Form.Sequence('axes', 'column labels of Euclidean axes',
                ('pc1', 'pc2', 'pc3')),
            Form.Float('radius',
                'remove duplicate points in this radius', '1e-8'),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    """
    @return: the format of the output
    """
    return FormOut.RTable()

def get_dup_indices(points, radius):
    """
    @param points: numpy points
    @param radius: do not allow points to be this close to each other
    """
    npoints = len(points)
    dups = set()
    for i in range(npoints):
        if i in dups:
            continue
        for j in range(i+1, npoints):
            if j in dups:
                continue
            d = np.linalg.norm(points[j] - points[i])
            if d < radius:
                dups.add(j)
    return dups

def get_response_content(fs):
    # read the table
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    Carbone.validate_headers(header_row)
    # get the numpy array of conformant points
    h_to_i = dict((h, i+1) for i, h in enumerate(header_row))
    axis_headers = fs.axes
    if not axis_headers:
        raise ValueError('no Euclidean axes were provided')
    axis_set = set(axis_headers)
    header_set = set(header_row)
    bad_axes = axis_set - header_set
    if bad_axes:
        raise ValueError('invalid axes: ' + ', '.join(bad_axes))
    axis_lists = []
    for h in axis_headers:
        index = h_to_i[h]
        try:
            axis_list = Carbone.get_numeric_column(data_rows, index)
        except Carbone.NumericError:
            msg_a = 'expected the axis column %s ' % h
            msg_b = 'to be numeric'
            raise ValueError(msg_a + msg_b)
        axis_lists.append(axis_list)
    points = np.array(zip(*axis_lists))
    # find the set of indices of duplicate points
    dup_indices = get_dup_indices(points, fs.radius)
    # get the data rows with duplicate indices removed
    new_rows = [row for i, row in enumerate(data_rows) if i not in dup_indices]
    # construct the new table
    out = StringIO()
    print >> out, '\t'.join(header_row)
    print >> out, '\n'.join('\t'.join(row) for row in new_rows)
    return out.getvalue()
