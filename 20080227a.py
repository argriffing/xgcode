"""Given N free energies, get an NxN rate matrix.
"""

from StringIO import StringIO
import math

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form
import iterutils

def get_form():
    """
    @return: the body of a form
    """
    # define the default energy string
    default_energies = [2, 4, 6, 8]
    default_energy_string = '\n'.join(str(E) for E in default_energies)
    # define the form objects
    form_objects = [
            Form.MultiLine('energies', 'ordered energies',
                default_energy_string)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the energies from the form data
    energies = []
    for line in iterutils.stripped_lines(StringIO(fs.energies)):
        try:
            energy = float(line)
        except ValueError, e:
            raise HandlingError('invalid energy: %s' % line)
        energies.append(energy)
    # create the row major matrix
    matrix = []
    # add the rates that are valid for off diagonal elements
    for row_energy in energies:
        row = []
        for column_energy in energies:
            rate = math.exp(-(column_energy - row_energy))
            row.append(rate)
        matrix.append(row)
    # correct the diagonal elements
    for i, row in enumerate(matrix):
        matrix[i][i] = -(sum(row) - 1)
    # create the string representing the rate matrix
    matrix_string = MatrixUtil.m_to_string(matrix)
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, matrix_string
