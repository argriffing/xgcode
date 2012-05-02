"""
Do some analysis of variance.

This uses the R software.
Also google for a pdf called Using R for an Analysis of Variance.
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

g_factor = 'cluster'

g_variable = 'temperature'


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_table),
            Form.SingleLine('factor', 'factor name', g_factor),
            Form.SingleLine('variable', 'variable name', g_variable)]
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
    if fs.variable not in header_row:
        msg = 'the variable name was not found as a column in the data table'
        raise ValueError(msg)
    if fs.factor not in header_row:
        msg = 'the factor name was not found as a column in the data table'
        raise ValueError(msg)
    return RUtil.run_with_table(fs.table, fs, get_script_content)

def get_script_content(fs, temp_table_name):
    """
    @param fs: something like a fieldstorage object
    @param temp_table_name: name of the temporary table file
    """
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'data <- data.frame(y=d$%s, group=factor(d$%s))' % (
                fs.variable, fs.factor),
            'myfit <- lm(y ~ group, data)',
            'anova(myfit)',
            'TukeyHSD(aov(y ~ group, data))']
    return '\n'.join(lines) + '\n'
