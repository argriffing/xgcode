"""
Plot endpoint-conditioned probabilities using weak evolution diffusion limits.

Plot endpoint-conditioned Wright-Fisher path state probabilities
using the diffusion limits of weak selection and mutation.
A simple two-allele model is used, but with mutation and selection.
The genome is assumed to have only a single site,
so the issue of recombination is moot.
If the diffusion approximation is good,
then it should not be too sensitive to the granularity parameter;
This insensitivity would mean that a large population could be
represented by a small state space without much loss of accuracy.
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
            Form.Float('pop', 'population size',
                '1e6', low_inclusive=1),
            Form.Integer('pop_gran',
                'diffusion approximation granularity (larger is more accurate)',
                40, low=2, high=200),
            Form.Integer('ngenerations', 'total number of generations',
                40, low=2, high=200),
            Form.Float('initial_freq', 'initial allele frequency',
                0.1, low_inclusive=0, high_inclusive=1),
            Form.Float('final_freq', 'final allele frequency',
                0.2, low_inclusive=0, high_inclusive=1),
            Form.Float('additive_selection', 'additive allele selection',
                '0.02'),
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
    # transform the arguments according to the diffusion approximation
    mutation_ab = (fs.pop * fs.mutation_ab) / fs.pop_gran
    mutation_ba = (fs.pop * fs.mutation_ba) / fs.pop_gran
    selection_ratio = 1 + (fs.pop * fs.additive_selection) / fs.pop_gran
    npop = fs.pop_gran
    ngenerations = fs.ngenerations
    nmutants_initial = int(fs.initial_freq * fs.pop_gran)
    nmutants_final = int(fs.final_freq * fs.pop_gran)
    # precompute some transition matrices
    P_drift_selection = pgmsinglesite.create_drift_selection_transition_matrix(
            npop, selection_ratio)
    MatrixUtil.assert_transition_matrix(P_drift_selection)
    P_mutation = pgmsinglesite.create_mutation_transition_matrix(
            npop, mutation_ab, mutation_ba)
    MatrixUtil.assert_transition_matrix(P_mutation)
    # define the R table headers
    headers = [
            'generation',
            'allele.frequency',
            'probability',
            'log.prob',
            ]
    # compute the transition matrix
    P = np.dot(P_drift_selection, P_mutation)
    # Compute the endpoint conditional probabilities for various states
    # along the unobserved path.
    nstates = npop + 1
    M = np.zeros((nstates, ngenerations))
    M[nmutants_initial, 0] = 1.0
    M[nmutants_final, ngenerations-1] = 1.0
    for i in range(ngenerations-2):
        A_exponent = i + 1
        B_exponent = ngenerations - 1 - A_exponent
        A = np.linalg.matrix_power(P, A_exponent)
        B = np.linalg.matrix_power(P, B_exponent)
        weights = np.zeros(nstates)
        for k in range(nstates):
            weights[k] = A[nmutants_initial, k] * B[k, nmutants_final]
        weights /= np.sum(weights)
        for k, p in enumerate(weights):
            M[k, i+1] = p
    arr = []
    for g in range(ngenerations):
        for k in range(nstates):
            p = M[k, g]
            if p:
                logp = math.log(p)
            else:
                logp = float('-inf')
            allele_frequency = k / float(npop)
            row = [g, allele_frequency, p, logp]
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

"""
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
"""

def get_ggplot():
    out = StringIO()
    print >> out, mk_call_str('require', '"reshape"')
    print >> out, mk_call_str('require', '"ggplot2"')
    print >> out, 'ggplot(data=my.table,'
    print >> out, mk_call_str(
            'aes', x='generation', y='allele.frequency', fill='log.prob')
    print >> out, ') + geom_tile() +',
    print >> out, mk_call_str(
            'scale_fill_continuous',
            breaks='c(-6, -5, -4, -3, -2, -1, 0)', limits='c(-6, 0)')
    return out.getvalue().rstrip()

