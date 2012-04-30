"""Compute the variation explained by a linear model.

This uses the R software.
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
                g_dependent_name)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # get the independent variable names
    indep = Util.get_stripped_lines(fs.independent.splitlines())
    dep = fs.dependent
    # get the r table
    rtable = RUtil.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    Carbone.validate_headers(header_row)
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
    return RUtil.run_with_table(fs.table, (indep, dep), get_script_content)

def get_script_content(data, temp_table_name):
    """
    @param data: the (indep, dep) data pair
    @param temp_table_name: name of the temporary table file
    """
    indep, dep = data
    symbolic_indep_sum = ' + '.join('d$' + x for x in indep)
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'myfit <- lm(d$%s ~ %s)' % (dep, symbolic_indep_sum),
            'summary(myfit)']
    return '\n'.join(lines) + '\n'
