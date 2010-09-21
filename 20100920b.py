"""Do some analysis of variance.

This uses the R software.
Also google for a pdf called Using R for an Analysis of Variance.
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

g_factor = 'cluster'

g_variable = 'temperature'


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_table),
            Form.SingleLine('factor', 'factor name', g_factor),
            Form.SingleLine('variable', 'variable name', g_variable),
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
    if fs.variable not in header_row:
        msg = 'the variable name was not found as a column in the data table'
        raise ValueError(msg)
    if fs.factor not in header_row:
        msg = 'the factor name was not found as a column in the data table'
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
            f_temp_table.name, fs.factor, fs.variable)
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

def get_script_content(temp_table_name, factor, variable):
    """
    @param temp_table_name: name of the temporary table file
    @param factor: a column name e.g. 'cluster'
    @param variable: a column name e.g. 'temperature'
    """
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'data <- data.frame(y=d$%s, group=factor(d$%s))' % (
                variable, factor),
            'myfit <- lm(y ~ group, data)',
            'anova(myfit)',
            'TukeyHSD(aov(y ~ group, data))']
    return '\n'.join(lines) + '\n'
