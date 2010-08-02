"""Given N free energies, get an NxN rate matrix.
"""

import math

from SnippetUtil import HandlingError
import MatrixUtil
import RateMatrix
import Form
import FormOut
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

def get_form_out():
    return FormOut.RateMatrix()

def get_response_content(fs):
    # read the energies from the form data
    energies = []
    for line in iterutils.stripped_lines(fs.energies.splitlines()):
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
    # return the rate matrix
    return MatrixUtil.m_to_string(matrix) + '\n'
