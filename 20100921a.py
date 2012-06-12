"""
Compute the Fisher exact test using R.
"""

from StringIO import StringIO
import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import RUtil
import Carbone
import const

g_tags = ['pca:compute']

g_table = const.read('20100920a')

g_var_a = 'cluster'
g_var_b = 'location'


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_table),
            Form.SingleLine('var_a', 'first variable', g_var_a),
            Form.SingleLine('var_b', 'second variable', g_var_b)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the r table
    rtable = RUtil.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    Carbone.validate_headers(header_row)
    # check requested variable names as column headers
    if fs.var_a not in header_row:
        raise ValueError('the first variable name is not column header')
    if fs.var_b not in header_row:
        raise ValueError('the second variable name is not column header')
    return RUtil.run_with_table(fs.table, fs, get_script_content)

def get_script_content(fs, temp_table_name):
    """
    @param fs: something like a fieldstorage object
    @param temp_table_name: name of the temporary table file
    """
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'fisher.test(factor(d$%s), factor(d$%s))' % (
                fs.var_a, fs.var_b)]
    return '\n'.join(lines) + '\n'
