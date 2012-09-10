"""
Check W-F mean and variance for the two compensatory substitution pathways.

Compare estimates obtained by forward simulation
to estimates obtained through numerical solutions of systems of equations
that approximate the forward simulation dynamics.
Two different compensatory substitution pathways are considered.
The model is of two alleles at two loci, so a four-haplotype model.
Mutation, recombination, and selection+drift occur between each generation
in that order.
Selection acts directly on the haploid population,
in a way that corresponds to multiplicative selection in diploids.
The two different compensatory substitution pathways
from the high-fitness state AB to the high-fitness state ab are distinguished
by whether the most recently fixed haplotype before first reaching state ab
is AB or is a low-fitness haplotype.
Haploid selection for this compensatory model
gives fitness 1 to AA and ab and fitness 1-s to Ab and aB.
"""

from StringIO import StringIO
import math

import numpy as np
from scipy import linalg

import Form
import FormOut
import RUtil
from RUtil import mk_call_str
import MatrixUtil
import wfengine
import wfcompens
import multinomstate
import wffwdcompens
import wfbckcompens
import combobreaker

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Integer('N_diploid', 'diploid population size N', 6,
                low=1, high=10),
            Form.Float('theta', 'mutation rate 4*N*mu', '1.0',
                low_inclusive=0.001),
            Form.Float('Nr', 'recombination rate N*r', '5.0',
                low_inclusive=0),
            Form.Float('Ns', 'selection N*s', '1.0',
                low_inclusive=0),
            Form.Integer(
                'nsamples_max',
                'max number of sampled paths', 500),
            ]

def get_form_out():
    return FormOut.Report()

def get_backward_info(N_diploid, theta, Nr, Ns):
    """
    Compute expectations and variances for the two substitution pathways.
    Here backward is somewhat of a misnomer; it is meant as a contrast
    to forward simulation.
    @param N_diploid: diploid population size
    @param theta: a mutation rate
    @param Nr: a recombination rate
    @param Ns: a selection value
    @return: (t1, v1), (t2, v2)
    """
    # set up the state space
    k = 4
    M = multinomstate.get_sorted_states(2*N_diploid, k)
    T = multinomstate.get_inverse_map(M)
    nstates = M.shape[0]
    lmcs = wfengine.get_lmcs(M)
    # compute rate matrices
    R_rate = wfcompens.create_recomb(M, T)
    M_rate = wfcompens.create_mutation(M, T)
    # compute a recombination probability matrix
    R_prob = linalg.expm(Nr * R_rate / float((2*N_diploid)**2))
    # compute the expected number of mutation events per generation
    mu = theta / 2
    # compute the mutation matrix
    # and the product of mutation and recombination.
    M_prob = linalg.expm(mu * M_rate / float(2*2*N_diploid))
    MR_prob = np.dot(M_prob, R_prob)
    # compute the selection coefficient
    s = Ns / float(N_diploid)
    lps = wfcompens.create_selection(s, M)
    S_prob = np.exp(wfengine.create_genic(lmcs, lps, M))
    P = np.dot(MR_prob, S_prob)
    #
    t1, v1 = wfbckcompens.get_type_1_info(P)
    t2, v2 = wfbckcompens.get_type_2_info(P)
    return (t1, v1), (t2, v2)

def sample_forward(N_hap, mu, r, fitnesses):
    """
    This uses the cython extension.
    @param N_hap: haploid population size or twice diploid population size
    @param mu: expected number of mutations per generation
    @param r: expected number of recombinations per generation
    @param fitnesses: relative haploid fitnesses of haplotypes
    @return: a list of meta states and dwell times
    """
    ntransitions = 0
    haplotype_probs = np.zeros(4, dtype=float)
    counts = np.array([N_hap, 0, 0, 0], dtype=int)
    state = np.empty(N_hap, dtype=int)
    #
    fixation_state = 0
    fixation_dwell = 0
    history = []
    #
    while True:
        # get the current population meta state
        current_fixation_state = -1
        for i in range(4):
            if counts[i] == N_hap:
                current_fixation_state = i
        # update the history list, fixation state, and dwell time
        if fixation_state == current_fixation_state:
            fixation_dwell += 1
        else:
            history.append((fixation_state, fixation_dwell))
            fixation_state = current_fixation_state
            fixation_dwell = 1
        # check if the compensatory substitution has occurred
        if fixation_state == 3:
            history.append((fixation_state, fixation_dwell))
            return history
        # expand the counts into the unary state
        wfcompens.expand_counts(state, counts)
        # do the mutation step
        nevents = np.random.poisson(mu)
        #if nevents > N_hap:
            #raise ValueError('sampling too many things without replacement')
        if nevents:
            individuals = np.random.randint(N_hap, size=nevents)
            loci = np.random.randint(2, size=nevents)
            wfcompens.multiple_mutation(state, counts, individuals, loci)
        # do the recombination step
        nevents = np.random.poisson(r)
        #if nevents*2 > N_hap:
            #raise ValueError('sampling too many things without replacement')
        if nevents:
            individuals = np.random.randint(N_hap, size=2*nevents)
            wfcompens.multiple_recombination(state, counts, individuals)
        # do the selection step
        wfcompens.reselection(haplotype_probs, fitnesses, counts)
        counts = np.random.multinomial(N_hap, haplotype_probs)
        # increment the number of observed generational transitions
        ntransitions += 1
    return ntransitions

class Accum:
    def __init__(self, N_hap, mu, r, s):
        self.N_hap = N_hap
        self.mu = mu
        self.r = r
        self.fitnesses = np.array([1, 1-s, 1-s, 1], dtype=float)
        self.type_1_times = []
        self.type_2_times = []
    def __call__(self):
        history = sample_forward(
                self.N_hap, self.mu, self.r, self.fitnesses)
        # Check history invariants.
        if history[0][0] != 0:
            raise ValueError('expected first state in the path to be fixed AB')
        if history[-1][0] != 3:
            raise ValueError('expected last state in the path to be fixed ab')
        # Compute the type 1 or type 2 wait time from the history.
        total_dwell = 0
        low_fitness_fixation = False
        for meta, dwell in reversed(history):
            total_dwell += dwell
            if meta in (1, 2):
                low_fitness_fixation = True
            if meta == 0:
                break
        if low_fitness_fixation:
            self.type_1_times.append(total_dwell)
        else:
            self.type_2_times.append(total_dwell)
        return False

def get_response_content(fs):
    #
    N_hap = fs.N_diploid * 2
    # TODO check s vs N_diploid so that selection is not too big
    # as to make fitnesses negative or zero
    #
    (t1, v1), (t2, v2) = get_backward_info(
            fs.N_diploid, fs.theta, fs.Nr, fs.Ns)
    #
    f = Accum(N_hap, fs.theta / 2, fs.Nr, fs.Ns / fs.N_diploid)
    #
    nseconds_max = 5.0
    info = combobreaker.run_callable(f, nseconds_max, fs.nsamples_max)
    #
    out = StringIO()
    print >> out, 'exact results:'
    print >> out
    print >> out, 'Type 1 compensatory substitution time expectation:'
    print >> out, t1
    print >> out
    print >> out, 'Type 1 compensatory substitution time variance:'
    print >> out, v1
    print >> out
    print >> out, 'Type 2 compensatory substitution time expectation:'
    print >> out, t2
    print >> out
    print >> out, 'Type 2 compensatory substitution time variance:'
    print >> out, v2
    print >> out
    print >> out
    print >> out, 'forward sampling results:'
    print >> out
    print >> out, 'Type 1 compensatory substitution time expectation:'
    print >> out, np.mean(f.type_1_times)
    print >> out
    print >> out, 'Type 1 compensatory substitution time variance:'
    print >> out, np.var(f.type_1_times)
    print >> out
    print >> out, 'Type 2 compensatory substitution time expectation:'
    print >> out, np.mean(f.type_2_times)
    print >> out
    print >> out, 'Type 2 compensatory substitution time variance:'
    print >> out, np.var(f.type_2_times)
    print >> out
    return out.getvalue()


def main():
    """
    Check some stuff that takes longer to compute.
    """
    pass

if __name__ == '__main__':
    main()

