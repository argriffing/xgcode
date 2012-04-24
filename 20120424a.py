"""
Sample some pairs of mutation and selection processes.

The idea is to make a csv file for viewing in ggobi,
maybe this will be useful or maybe not.
"""

from StringIO import StringIO
import random
import math
from math import log

import numpy as np

import Form
import FormOut
import mrate
import ctmcmi
import divtime
import combobreaker

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('nsamples', 'sample at most this many data points',
                1000, low=1)]
    return form_objects

def get_form_out():
    return FormOut.Csv('mut-sel-samples')

class Accumulator:
    def __init__(self):
        self.rows = []
    def __call__(self):
        self.rows.append(sample_row())

def sample_row():
    n = 4
    # sample the exchangeability
    S = np.zeros((n, n))
    S[1,0] = random.expovariate(1)
    S[2,0] = random.expovariate(1)
    S[2,1] = random.expovariate(1)
    S[3,0] = random.expovariate(1)
    S[3,1] = random.expovariate(1)
    S[3,2] = random.expovariate(1)
    # sample the mutation stationary distribution
    mdistn = np.array([random.expovariate(1) for i in range(n)])
    mdistn /= np.sum(mdistn)
    # sample the mutation selection balance stationary distribution
    bdistn = np.array([random.expovariate(1) for i in range(n)])
    bdistn /= np.sum(bdistn)
    # sample the time
    t = random.expovariate(1)
    # sample the info type
    infotype = random.choice(('infotype.mi', 'infotype.fi'))
    # Compute some intermediate variables
    # from which the summary statistics and the label are computed.
    S = S + S.T
    M = S * mdistn
    M -= np.diag(np.sum(M, axis=1))
    R = mrate.to_gtr_halpern_bruno(M, bdistn)
    shannon_ent_mut = -sum(p*log(p) for p in mdistn)
    shannon_ent_bal = -sum(p*log(p) for p in bdistn)
    logical_ent_mut = 1.0 - sum(p*p for p in mdistn)
    logical_ent_bal = 1.0 - sum(p*p for p in bdistn)
    expected_rate_mut = mrate.Q_to_expected_rate(M)
    expected_rate_bal = mrate.Q_to_expected_rate(R)
    spectral_rate_mut = 1 / mrate.R_to_relaxation_time(M)
    spectral_rate_bal = 1 / mrate.R_to_relaxation_time(R)
    mi_mut = ctmcmi.get_mutual_information(M, t)
    mi_bal = ctmcmi.get_mutual_information(R, t)
    fi_mut = divtime.get_fisher_information(M, t)
    fi_bal = divtime.get_fisher_information(R, t)
    # compute the summary statistics
    summary_entries = [
            shannon_ent_bal - shannon_ent_mut,
            logical_ent_bal - logical_ent_mut,
            log(shannon_ent_bal) - log(shannon_ent_mut),
            log(logical_ent_bal) - log(logical_ent_mut),
            expected_rate_bal - expected_rate_mut,
            spectral_rate_bal - spectral_rate_mut,
            log(expected_rate_bal) - log(expected_rate_mut),
            log(spectral_rate_bal) - log(spectral_rate_mut),
            mi_bal - mi_mut,
            fi_bal - fi_mut,
            math.log(mi_bal) - math.log(mi_mut),
            math.log(fi_bal) - math.log(fi_mut),
            ]
    # get the definition entries
    definition_entries = [
            S[1,0], S[2,0], S[2,1], S[3,0], S[3,1], S[3,2],
            mdistn[0], mdistn[1], mdistn[2], mdistn[3],
            bdistn[0], bdistn[1], bdistn[2], bdistn[3],
            infotype,
            t,
            ]
    # define the label
    if infotype == 'infotype.mi' and mi_mut > mi_bal:
        label = 'mut.is.better'
    elif infotype == 'infotype.mi' and mi_mut < mi_bal:
        label = 'bal.is.better'
    elif infotype == 'infotype.fi' and fi_mut > fi_bal:
        label = 'mut.is.better'
    elif infotype == 'infotype.fi' and fi_mut < fi_bal:
        label = 'bal.is.better'
    else:
        label = 'indistinguishable'
    # return the row
    return definition_entries + summary_entries + [label]

def get_response_content(fs):
    headers = [
            # define the point
            'exch21', 'exch31', 'exch32', 'exch41', 'exch42', 'exch43',
            'mutprob1', 'mutprob2', 'mutprob3', 'mutprob4',
            'balprob1', 'balprob2', 'balprob3', 'balprob4',
            'infotype',
            't',
            # annotate with some summary statistics
            'shannon.ent.diff',
            'logical.ent.diff',
            'log.shannon.ent.ratio',
            'log.logical.ent.ratio',
            'expected.rate.diff',
            'spectral.rate.diff',
            'log.expected.rate.ratio',
            'log.spectral.rate.ratio',
            'mi.diff',
            'fi.diff',
            'log.mi.ratio',
            'log.fi.ratio',
            # define the label
            'label',
            ]
    accum = Accumulator()
    combobreaker.run_callable(accum, nseconds=5, niterations=fs.nsamples)
    # return csv contents
    out = StringIO()
    print >> out, ', '.join(headers)
    for row in accum.rows:
        print >> out, ', '.join(str(x) for x in row)
    return out.getvalue()

