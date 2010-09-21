"""Compute the Fisher exact test using R.
"""

from StringIO import StringIO
import os
import tempfile

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
            Form.SingleLine('var_b', 'second variable', g_var_b),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the r table
    rtable = Carbone.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    # Do a more stringent check of the column headers.
    for h in header_row:
        if not Carbone.is_valid_header(h):
            msg = 'invalid column header: %s' % h
            raise ValueError(msg)
    # check requested variable names as column headers
    if fs.var_a not in header_row:
        msg = 'the first variable name is not column header'
        raise ValueError(msg)
    if fs.var_b not in header_row:
        msg = 'the second variable name is not column header'
        raise ValueError(msg)
    # define the temp table content
    temp_table_content = fs.table
    # Create a temporary data table file for R.
    f_temp_table = tempfile.NamedTemporaryFile(delete=False)
    f_temp_table.write(temp_table_content)
    f_temp_table.close()
    # Create a temporary pathname for the plot created by R.
    temp_plot_name = Util.get_tmp_filename()
    # Create a temporary R script file.
    f_temp_script = tempfile.NamedTemporaryFile(delete=False)
    script_content = get_script_content(
            f_temp_table.name, fs.var_a, fs.var_b)
    f_temp_script.write(script_content)
    f_temp_script.close()
    # Call R.
    retcode, r_out, r_err = RUtil.run(f_temp_script.name)
    if retcode:
        raise ValueError('R error:\n' + r_err)
    # Delete the temporary data table file.
    os.unlink(f_temp_table.name)
    # Delete the temporary script file.
    os.unlink(f_temp_script.name)
    # Return the R stderr as a string.
    return r_err

def get_script_content(temp_table_name, var_a, var_b):
    """
    @param temp_table_name: name of the temporary table file
    @param var_a: the first column name
    @param var_b: the second column name
    """
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'fisher.test(factor(d$%s), factor(d$%s))' % (
                var_a, var_b)]
    return '\n'.join(lines) + '\n'
