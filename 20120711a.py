"""
Sample some endpoint conditioned mutant allele frequency histories.

Everything about this model is discrete and finite;
diffusion limits are not considered.
The selection parameter defines the ratio of the probability that a
given mutant will be Wright-Fisher selected for the next generation
to the probability that a given wild-type will be selected.
"""

from StringIO import StringIO

import numpy as np

import Form
import FormOut
import PathSampler
import MatrixUtil
import pgmsinglesite
import RUtil
from RUtil import mk_call_str


def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('npop', 'population size',
                100, low=1, high=200),
            Form.Integer('ngenerations', 'total number of generations',
                100, low=2, high=200),
            Form.Integer('nmutants_initial', 'initial number of mutants',
                10, low=0),
            Form.Integer('nmutants_final', 'final number of mutants',
                20, low=0),
            Form.Float('selection_ratio', 'selection probability ratio',
                '2.0', low_exclusive=0),
            Form.Float('mutation_ab',
                'wild-type to mutant transition probability',
                '0.01', low_inclusive=0, high_inclusive=1),
            Form.Float('mutation_ba',
                'mutant to wild-type transition probability',
                '0.01', low_inclusive=0, high_inclusive=1),
            Form.ImageFormat(),
            ]

def get_form_out():
    return FormOut.Image('plot')

def get_response_content(fs):
    # precompute some transition matrices
    P_drift_selection = pgmsinglesite.create_drift_selection_transition_matrix(
            fs.npop, fs.selection_ratio)
    MatrixUtil.assert_transition_matrix(P_drift_selection)
    P_mutation = pgmsinglesite.create_mutation_transition_matrix(
            fs.npop, fs.mutation_ab, fs.mutation_ba)
    MatrixUtil.assert_transition_matrix(P_mutation)
    # define the R table headers
    headers = ['generation', 'number.of.mutants']
    # compute the path samples
    P = np.dot(P_drift_selection, P_mutation)
    mypath = PathSampler.sample_endpoint_conditioned_path(
            fs.nmutants_initial, fs.nmutants_final, fs.ngenerations, P)
    arr = [[i, nmutants] for i, nmutants in enumerate(mypath)]
    # create the R table string and scripts
    # get the R table
    table_string = RUtil.get_table_string(arr, headers)
    # get the R script
    script = get_ggplot()
    # create the R plot image
    device_name = Form.g_imageformat_to_r_function[fs.imageformat]
    retcode, r_out, r_err, image_data = RUtil.run_plotter(
            table_string, script, device_name)
    if retcode:
        raise RUtil.RError(r_err)
    return image_data

def get_ggplot():
    out = StringIO()
    print >> out, mk_call_str('require', '"reshape"')
    print >> out, mk_call_str('require', '"ggplot2"')
    print >> out, 'my.table.long <-',
    print >> out, mk_call_str('melt', 'my.table', id='"generation"')
    print >> out, 'ggplot(data=my.table.long,'
    print >> out, mk_call_str(
            'aes', x='generation', y='value', colour='variable')
    print >> out, ') + geom_line()',
    return out.getvalue()

