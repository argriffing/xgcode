"""Given a rate matrix, check the detailed balance equations for reversibility.
"""

import StringIO
import math

import numpy

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the default rate matrix
    dictionary_rate_matrix = RateMatrix.get_sample_codon_rate_matrix()
    labels = list(sorted(set(a for a, b in dictionary_rate_matrix)))
    R = numpy.array(MatrixUtil.dict_to_row_major(dictionary_rate_matrix, labels, labels))
    return [Form.Matrix('matrix', 'rate matrix', R, MatrixUtil.assert_rate_matrix)]

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the matrix from the form data
    R = fs.matrix
    # get the stationary distribution of the rate matrix
    try:
        stationary_distribution = RateMatrix.get_stationary_distribution(R.tolist())
    except RateMatrix.RateMatrixError, e:
        raise HandlingError('error calculating the stationary distribution: ' + str(e))
    # for each pair of entries, check the detailed balance equation
    table_rows = []
    for i, pi_i in enumerate(stationary_distribution):
        for j, pi_j in enumerate(stationary_distribution):
            r_ij = R[i][j]
            r_ji = R[j][i]
            if pi_i*r_ij != pi_j*r_ji:
                row = []
                row.append(abs(math.log(pi_i * r_ij) - math.log(pi_j * r_ji)))
                row.extend([pi_i, pi_j, r_ij, r_ji])
                table_rows.append(row)
    # write some stuff
    out = StringIO.StringIO()
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
        response_headers = [('Content-Type', 'text/html')]
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
        response_headers = [('Content-Type', 'text/plain')]
        print >> out, 'All detailed balance equations are satisfied for this rate matrix.'
    # return the response
    return response_headers, out.getvalue().strip()
