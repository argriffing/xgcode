"""
Plot a clustering of columns of an R table using squared correlation.

This uses the R software.
Correlation is rounded to two decimal places,
and if correlation is not defined then the value .01 is used.
Correlation is computed using the Kendall method.
Clustering is done using the hclust function in R and
with distances equal to one minus the squared correlation.
"""

import os

from SnippetUtil import HandlingError
import Form
import FormOut
import RUtil
import Util

g_tags = ['pca:plot']

g_rows = [
        ('IC31', 'IC32', 'IC33'),
        (0, 2, 0, 0),
        (1, 0, 2, 0),
        (2, 0, 0, 2),
        (3, 2, 0, 0),
        (4, 0, 2, 0),
        (5, 0, 0, 2),
        (6, 2, 0, 0),
        (7, 0, 2, 0),
        (8, 0, 0, 1),
        (9, 0, 0, 1),
        (10, 1, 1, 0),
        (11, 1, 1, 0),
        (12, 0, 0, 2),
        (13, 2, 2, 0),
        (14, 0, 0, 1),
        (15, 0, 0, 1),
        (16, 2, 2, 0),
        (17, 0, 0, 2)]

g_lines = ['\t'.join(str(x) for x in row) for row in g_rows]

g_script_body = """
RR <- round(cor(my.table)^2, 2)
RR[is.na(RR) == TRUE] <- 0.01
h <- hclust(as.dist(1 - RR))
plot(h)
"""

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'table', '\n'.join(g_lines)),
            Form.ImageFormat()]
    return form_objects

def get_form_out():
    return FormOut.Image('hclust')

def get_response_content(fs):
    # check the r table
    RUtil.RTable(fs.table.splitlines())
    # make the plot
    device = Form.g_imageformat_to_r_function[fs.imageformat]
    image_data = RUtil.run_plotter_concise(
            fs.table, g_script_body, device)
    return image_data
