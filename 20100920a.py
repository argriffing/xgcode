"""Compute the variation explained by a linear model.

This uses the R software.
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

g_table = const.read('20100709a')

g_independent_names = ['temperature', 'precipitation']

g_dependent_name = 'pc1'


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_table),
            Form.MultiLine('independent', 'names of independent variables',
                '\n'.join(g_independent_names)),
            Form.SingleLine('dependent', 'name of the dependent variable',
                g_dependent_name),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the independent variable names
    indep = Util.get_stripped_lines(fs.independent.splitlines())
    dep = fs.dependent
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
    bad_indep_names = set(indep) - set(header_row)
    if bad_indep_names:
        msg_a = 'these requested independent variable names '
        msg_b = 'were not found as columns in the data table: '
        msg_c = str(bad_indep_names)
        raise ValueError(msg_a + msg_b + msg_c)
    if dep not in header_row:
        msg_a = 'the dependent variable name '
        msg_b = 'was not found as a column in the data table'
        raise ValueError(msg_a + msg_b)
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
    script_content = get_script_content(f_temp_table.name, indep, dep)
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

def get_script_content(temp_table_name, indep, dep):
    """
    @param temp_table_name: name of the temporary table file
    @param indep: list of the independent variable names
    @param dep: the dependent variable name
    """
    symbolic_indep_sum = ' + '.join('d$' + x for x in indep)
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'myfit <- lm(d$%s ~ %s)' % (dep, symbolic_indep_sum),
            'summary(myfit)',
            'cat("ohai")']
    return '\n'.join(lines) + '\n'
