"""
Plot endpoint-conditioned Wright-Fisher path state probabilities.

Although populations can be large-ish,
the size of genome is always a single position with two alleles.
Everything about this model is discrete and finite;
diffusion limits are not considered.
The selection parameter defines the ratio of the probability that a
given mutant will be Wright-Fisher selected for the next generation
to the probability that a given wild-type will be selected.
"""

from StringIO import StringIO
import math

import numpy as np

import Form
import FormOut
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
                40, low=1, high=200),
            Form.Integer('ngenerations', 'total number of generations',
                40, low=2, high=200),
            Form.Integer('nmutants_initial', 'initial number of mutants',
                5, low=0),
            Form.Integer('nmutants_final', 'final number of mutants',
                10, low=0),
            Form.Float('selection_ratio', 'selection probability ratio',
                '1.0', low_exclusive=0),
            Form.Float('mutation_ab',
                'wild-type to mutant transition probability',
                '0.02', low_inclusive=0, high_inclusive=1),
            Form.Float('mutation_ba',
                'mutant to wild-type transition probability',
                '0.02', low_inclusive=0, high_inclusive=1),
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
    headers = [
            'generation',
            'number.of.mutants',
            'probability',
            'log.prob',
            ]
    # compute the transition matrix
    P = np.dot(P_drift_selection, P_mutation)
    # Compute the endpoint conditional probabilities for various states
    # along the unobserved path.
    nstates = fs.npop + 1
    M = np.zeros((nstates, fs.ngenerations))
    M[fs.nmutants_initial, 0] = 1.0
    M[fs.nmutants_final, fs.ngenerations-1] = 1.0
    for i in range(fs.ngenerations-2):
        A_exponent = i + 1
        B_exponent = fs.ngenerations - 1 - A_exponent
        A = np.linalg.matrix_power(P, A_exponent)
        B = np.linalg.matrix_power(P, B_exponent)
        weights = np.zeros(nstates)
        for k in range(nstates):
            weights[k] = A[fs.nmutants_initial, k] * B[k, fs.nmutants_final]
        weights /= np.sum(weights)
        for k, p in enumerate(weights):
            M[k, i+1] = p
    arr = []
    for g in range(fs.ngenerations):
        for k in range(nstates):
            p = M[k, g]
            if p:
                logp = math.log(p)
            else:
                logp = float('-inf')
            row = [g, k, p, logp]
            arr.append(row)
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

def get_presets():
    return [
            Form.Preset(
                'bimodal path distribution',
                {
                    'npop' : '40',
                    'ngenerations' : '40',
                    'nmutants_initial' : '5',
                    'nmutants_final' : '10',
                    'selection_ratio' : '2.0',
                    'mutation_ab' : '0.02',
                    'mutation_ba' : '0.02',
                    }),
            Form.Preset(
                'unimodal high',
                {
                    'npop' : '40',
                    'ngenerations' : '40',
                    'nmutants_initial' : '5',
                    'nmutants_final' : '10',
                    'selection_ratio' : '2.0',
                    'mutation_ab' : '0.1',
                    'mutation_ba' : '0.1',
                    }),
            Form.Preset(
                'unimodal low',
                {
                    'npop' : '40',
                    'ngenerations' : '40',
                    'nmutants_initial' : '5',
                    'nmutants_final' : '10',
                    'selection_ratio' : '2.0',
                    'mutation_ab' : '0.01',
                    'mutation_ba' : '0.01',
                    }),
                ]

def get_ggplot():
    out = StringIO()
    print >> out, mk_call_str('require', '"reshape"')
    print >> out, mk_call_str('require', '"ggplot2"')
    #print >> out, 'my.table.long <-',
    #print >> out, mk_call_str('melt', 'my.table', id='"generation"')
    #print >> out, 'ggplot(data=my.table.long,'
    print >> out, 'ggplot(data=my.table,'
    print >> out, mk_call_str(
            'aes', x='generation', y='number.of.mutants', fill='log.prob')
    print >> out, ') + geom_tile()',
    #print >> out, '+ scale_fill_continuous(breaks=c(0,0.001,0.01,0.1,1.0))',
    print >> out, '+ scale_fill_continuous(',
    print >> out, 'breaks=c(-6, -5, -4, -3, -2, -1, 0),',
    print >> out, 'limits=c(-6, 0)',
    print >> out, ')',
    return out.getvalue()

