"""Given a rate matrix, check the detailed balance equations for reversibility.
"""

from StringIO import StringIO
import math

import numpy as np

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the default rate matrix
    dictionary_rate_matrix = RateMatrix.get_sample_codon_rate_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    R = MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels)
    R = np.array(R)
    form_objects = [
            Form.Matrix('matrix', 'rate matrix',
                R, MatrixUtil.assert_rate_matrix)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the matrix from the form data
    R = fs.matrix
    # get the stationary distribution of the rate matrix
    try:
        v = RateMatrix.get_stationary_distribution(R.tolist())
    except RateMatrix.RateMatrixError as e:
        msg = 'error calculating the stationary distribution: ' + str(e)
        raise HandlingError(msg)
    # for each pair of entries, check the detailed balance equation
    table_rows = []
    for i, pi_i in enumerate(v):
        for j, pi_j in enumerate(v):
            r_ij = R[i][j]
            r_ji = R[j][i]
            if pi_i*r_ij != pi_j*r_ji:
                row = []
                row.append(abs(math.log(pi_i * r_ij) - math.log(pi_j * r_ji)))
                row.extend([pi_i, pi_j, r_ij, r_ji])
                table_rows.append(row)
    # write some stuff
    out = StringIO()
    if table_rows:
        # get the detailed balance html rows
        detailed_balance_rows = []
        for row in reversed(list(sorted(table_rows))):
            detailed_balance_rows.append(''.join('<td>' + str(value) + '</td>' for value in row))
        # get the header row
        header_entries = []
        header_entries.append('abs(log(&pi;<sub>i</sub>r<sub>ij</sub>)-log(&pi;<sub>j</sub>r<sub>ji</sub>))')
        header_entries.append('&pi;<sub>i</sub>')
        header_entries.append('&pi;<sub>j</sub>')
        header_entries.append('r<sub>ij</sub>')
        header_entries.append('r<sub>ji</sub>')
        header_row = ''.join('<th>%s</th>' % entry for entry in header_entries)
        # show detailed balance equation results
        print >> out, '<p>'
        print >> out, 'This table shows each state pair for which the detailed balance equation is not satisfied exactly.'
        print >> out, '</p>'
        print >> out, '<html>'
        print >> out, '<body>'
        print >> out, '<table>'
        print >> out, '<tr>' + header_row + '</tr>'
        for row in detailed_balance_rows:
            print >> out, '<tr>' + row + '</tr>'
        print >> out, '</table>'
        print >> out, '</body>'
        print >> out, '</html>'
    else:
        print >> out, '<html><body>'
        print >> out, 'All detailed balance equations are satisfied for this rate matrix.'
        print >> out, '</body></html>'
    # return the response
    return out.getvalue()
