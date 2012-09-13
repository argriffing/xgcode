"""
Check W-F ancestral vs. population process compensatory pathway preference.

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
import wffwdkline
import combobreaker
import iterutils

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('theta', 'mutation rate 4*N*mu', '1.0',
                low_inclusive=0.001),
            Form.Float('Ns', 'selection N*s', '1.0',
                low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report()

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

class AccumAncestral():
    def __init__(self, N_hap, mu, s):
        self.N_hap = N_hap
        self.mu = mu
        self.fitnesses = np.array([1, 1-s, 1-s, 1], dtype=float)
        self.generation = 0
        # define the original background state
        parent = None
        state = wffwdkline.g_AB
        m = wffwdkline.Mutation(self.generation, parent, state)
        m.count = N_hap
        #
        self.fitness_history = [self.fitnesses[state]]
        self.background = m
        self.mutations = [m]
    def __call__(self):
        """
        Each call is a generation.
        """
        self.mutations = wffwdkline.evolve(
                self.N_hap,
                self.mu,
                self.generation,
                self.mutations,
                self.fitnesses,
                self.fitness_history,
                )
        self.generation += 1
        return False

def extract_ancestral_lineage(background):
    """
    @param background: the original mutational path
    @return: the sequence of mutational paths to the mrca
    """
    m = background
    ancestral_lineage = []
    while True:
        if len(m.child_id_map) > 1 or m.count:
            return ancestral_lineage
        else:
            child = m.child_id_map.values()[0]
            ancestral_lineage.append(child)
            m = child
    raise ValueError('no ancestral line was found')

def extract_type_counts(ancestral_lineage):
    """
    @return: type 1 substitution counts, type 2 substitution counts
    """
    high_fitness_states = (wffwdkline.g_AB, wffwdkline.g_ab)
    low_fitness_states = (wffwdkline.g_Ab, wffwdkline.g_aB)
    nt1 = 0
    nt2 = 0
    prev_high_fitness_state = None
    prev_state = None
    for m in ancestral_lineage:
        # check for a compensatory subsitution on the ancestral lineage
        if m.state == wffwdkline.g_ab:
            if prev_high_fitness_state == wffwdkline.g_AB:
                if prev_state in low_fitness_states:
                    nt1 += 1
                elif prev_state == wffwdkline.g_AB:
                    nt2 += 1
        # record some historical context of the ancestral lineage
        if m.state in high_fitness_states:
            prev_high_fitness_state = m.state
        prev_state = m.state
    # return the counts of the compensatory substitution types
    return nt1, nt2

def get_ancestral_excess_fitness(ancestral_lineage, fitnesses, fitness_history):
    excess_total = 0
    for ma, mb in iterutils.pairwise(ancestral_lineage):
        for generation in range(ma.generation, mb.generation):
            excess_total += fitnesses[ma.state] - fitness_history[generation]
    #
    ngenerations_total = 0
    ngenerations_total += ancestral_lineage[-1].generation
    ngenerations_total -= ancestral_lineage[0].generation
    #
    return excess_total / ngenerations_total

def get_response_content(fs):
    # hardcode some values
    N_diploid = 100
    N_hap = N_diploid * 2
    nseconds = 3
    Nr = 0
    #
    out = StringIO()
    print >> out, 'theta:', fs.theta
    print >> out, 'Ns:', fs.Ns
    print >> out, '2N:', N_hap
    print >> out
    #
    f = Accum(N_hap, fs.theta / 2, Nr, fs.Ns / N_diploid)
    info = combobreaker.run_callable(f, nseconds=nseconds)
    nt1 = len(f.type_1_times)
    nt2 = len(f.type_2_times)
    print >> out, 'fixed state in the population:'
    print >> out, 'number of type 1 substitutions:', nt1
    print >> out, 'number of type 2 substitutions:', nt2
    pt1 = nt1 / float(nt1 + nt2)
    pt2 = nt2 / float(nt1 + nt2)
    print >> out, 'proportion of type 1 substitutions:', pt1
    print >> out, 'proportion of type 2 substitutions:', pt2
    print >> out
    f = AccumAncestral(N_hap, fs.theta / 2, fs.Ns / N_diploid)
    info = combobreaker.run_callable(f, nseconds=nseconds)
    ancestral_lineage = extract_ancestral_lineage(f.background)
    """
    for m in ancestral_lineage:
        print m.generation, m.state
    print
    """
    nt1, nt2 = extract_type_counts(ancestral_lineage)
    print >> out, 'state of the ancestral lineage:'
    print >> out, 'number of type 1 substitutions:', nt1
    print >> out, 'number of type 2 substitutions:', nt2
    pt1 = nt1 / float(nt1 + nt2)
    pt2 = nt2 / float(nt1 + nt2)
    print >> out, 'proportion of type 1 substitutions:', pt1
    print >> out, 'proportion of type 2 substitutions:', pt2
    print >> out
    s = fs.Ns / N_diploid
    fitnesses = np.array([1, 1-s, 1-s, 1], dtype=float)
    print >> out, 'excess fitness of fitter allele:', s
    print >> out, 'excess fitness of ancestral lineage:',
    print >> out, get_ancestral_excess_fitness(
            ancestral_lineage, fitnesses, f.fitness_history)
    print >> out
    return out.getvalue()


def main():
    """
    Check some stuff that takes longer to compute.
    """
    pass

if __name__ == '__main__':
    main()

