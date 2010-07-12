"""Combine R tables.
"""

from StringIO import StringIO
import os

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone
import iterutils


g_default_info_rows = [
        ['temperature', 'precipitation', 'otu'],
        ['1', '15', '600', 'IC31'],
        ['2', '30', '700', 'IC32'],
        ['3', '45', '800', 'IC33']]

g_default_info_lines = ['\t'.join(row) for row in g_default_info_rows]


g_default_pca_rows = [
        ['otu', 'pc1', 'pc2', 'pc3'],
        ['1', 'IC32', '0.563951577482', '0', '-0.426566077245'],
        ['2', 'IC31', '1.48488064713', '0', '0.270013409975'],
        ['3', 'IC33', '-1.0244161123', '-0.942809041582', '0.0782763336347']]

g_default_pca_lines = ['\t'.join(row) for row in g_default_pca_rows]


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table_a', 'first R table',
                '\n'.join(g_default_info_lines)),
            Form.MultiLine('table_b', 'second R table',
                '\n'.join(g_default_pca_lines)),
            Form.SingleLine('join_header', 'combine using this column', 'otu'),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable('out.table', [])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs, fs.table_a.splitlines(), fs.table_b.splitlines())
    filename = 'fungus.table'
    response_headers = [('Content-Type', 'text/plain')] 
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers.append(('Content-Disposition', disposition)) 
    return response_headers, text

def process(args, raw_a_lines, raw_b_lines):
    a_table = Carbone.RTable(raw_a_lines)
    b_table = Carbone.RTable(raw_b_lines)
    if args.join_header not in a_table.headers:
        msg = 'the first table does not have the requested column'
        raise ValueError(msg)
    if args.join_header not in b_table.headers:
        msg = 'the second table does not have the requested column'
        raise ValueError(msg)
    concat_headers = a_table.headers + b_table.headers
    out_headers = list(iterutils.unique_everseen(concat_headers))
    nunique_headers = len(out_headers)
    if len(a_table.headers) + len(b_table.headers) != nunique_headers + 1:
        msg = 'the tables should share only the requested column'
        raise ValueError(msg)
    # get the column index for each table
    a_index = a_table.header_to_column_index(args.join_header)
    b_index = b_table.header_to_column_index(args.join_header)
    # get the join column for each table
    a_column = a_table.header_to_primary_column(args.join_header)
    b_column = b_table.header_to_primary_column(args.join_header)
    # get the set of join elements common to both tables
    common_join_element_set = set(a_column) & set(b_column)
    # for the second table get the map from join elements to row indices
    b_j_to_i = dict((j, i) for i, j in enumerate(b_column))
    # create the output table without the R row labels
    out_data = []
    out_r_label = 1
    for row in a_table.data:
        a_join_element = row[a_index]
        if a_join_element not in common_join_element_set:
            continue
        a_out_row = row[1:]
        b_row_index = b_j_to_i[a_join_element]
        b_row = b_table.data[b_row_index]
        b_out_row = [x for i, x in enumerate(b_row) if i != b_index][1:]
        out_row = [out_r_label] + a_out_row + b_out_row
        out_data.append(out_row)
        out_r_label += 1
    # write the R table
    out = StringIO()
    print >> out, '\t'.join(out_headers)
    for row in out_data:
        print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue()
