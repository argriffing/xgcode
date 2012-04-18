"""Plot a clustering of columns of an R table using squared correlation.

This uses the R software.
Correlation is rounded to two decimal places,
and if correlation is not defined then the value .01 is used.
Correlation is computed using the Kendall method.
Clustering is done using the hclust function in R and
with distances equal to one minus the squared correlation.
"""

import tempfile
import os

from SnippetUtil import HandlingError
import Form
import FormOut
import hud
import Phylip
import RUtil
import Util
import Carbone

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
    # get the r table
    rtable = RUtil.RTable(fs.table.splitlines())
    header_row = rtable.headers
    data_rows = rtable.data
    # Do a more stringent check of the column headers.
    for h in header_row:
        if not Carbone.is_valid_header(h):
            msg = 'invalid column header: %s' % h
            raise ValueError(msg)
    # define the temp table content and the R image format function
    temp_table_content = fs.table
    image_function_name = Form.g_imageformat_to_r_function[fs.imageformat]
    # Create a temporary data table file for R.
    f_temp_table = tempfile.NamedTemporaryFile(delete=False)
    f_temp_table.write(temp_table_content)
    f_temp_table.close()
    # Create a temporary pathname for the plot created by R.
    temp_plot_name = Util.get_tmp_filename()
    # Create a temporary R script file.
    f_temp_script = tempfile.NamedTemporaryFile(delete=False)
    script_content = get_script_content(
            temp_plot_name, f_temp_table.name, image_function_name)
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
    # Read the image file.
    try:
        with open(temp_plot_name, 'rb') as fin:
            image_data = fin.read()
    except IOError as e:
        raise HandlingError('the R call seems to not have created the plot')
    # Delete the temporary image file.
    os.unlink(temp_plot_name)
    # Return the image data as a string.
    return image_data

def get_script_content(temp_plot_name, temp_table_name, image_function_name):
    lines = [
            'd <- read.table("%s")' % temp_table_name,
            'RR <- round(cor(d)^2, 2)',
            'RR[is.na(RR) == TRUE] <- 0.01',
            'h <- hclust(as.dist(1 - RR))',
            '%s("%s")' % (image_function_name, temp_plot_name),
            'plot(h)',
            'dev.off()']
    return '\n'.join(lines) + '\n'
