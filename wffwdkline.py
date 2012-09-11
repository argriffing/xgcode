"""
Implement a slow klineage-like Wright-Fisher forward simulation.

This module is loosely inspired by klineage,
an stl+boost+gsl software by Reed Cartwright.
For generality, a few details should be left to the caller of the module.
The mutation transition matrix should be provided,
as should the fitnesses.
The first intended application of this module is to compensatory evolution.
This is an evolutionary model with four microstates whose mutation graph
is like a square and where two opposite corners have high fitness
and the complementary opposite corners have low fitness.
So all mutations occur between high and low fitness states.
The mutation happens before the selection+drift.
Mutations are tracked by descent rather than by state.
It may be more natural to have mutation objects,
but I will try using a less object oriented approach.
"""

import unittest
import random

import numpy as np


g_AB = 0
g_Ab = 1
g_aB = 2
g_ab = 3

g_mutation = {
        (g_AB, 0) : g_aB,
        (g_AB, 1) : g_Ab,
        (g_Ab, 0) : g_ab,
        (g_Ab, 1) : g_AB,
        (g_aB, 0) : g_AB,
        (g_aB, 1) : g_ab,
        (g_ab, 0) : g_Ab,
        (g_ab, 1) : g_aB,
        }

class Mutation:
    def __init__(self, generation, parent, state):
        """
        @param generation: the generation at which the mutation arose
        @param parent: a link to the parent mutation object
        @param state: the haplotype state
        """
        self.generation = generation
        self.parent = parent
        self.state = state
        self.count = 1
        self.child_id_map = {}
    def is_dead(self):
        return not self.count and not self.child_id_map

def mutate(mu, generation, mutations):
    """
    Modify mutations in-place.
    @param mu: expected total mutations per generation
    @param generation: the current generation
    @param mutations: a sequence of mutation objects to be updated
    @return: None
    """
    # sample the number of mutations that occur at this generation
    nevents = np.random.poisson(mu)
    if not nevents:
        return
    #
    expanded = []
    for m in mutations:
        expanded.extend([m]*m.count)
    #
    replaced_mutations = random.sample(expanded, nevents)
    # create the new mutations
    w = np.random.randint(2, size=nevents)
    for i, m in enumerate(replaced_mutations):
        m = replaced_mutations[i]
        new_mutation = Mutation(generation, m, g_mutation[m.state, w[i]])
        m.count -= 1
        m.child_id_map[id(new_mutation)] = new_mutation
        mutations.append(new_mutation)
    return

def reselect(N_hap, mutations, fitnesses):
    """
    @param N_hap: total haplotype population size
    @param mutations: a sequence of mutation objects to be updated
    @param fitnesses: fitness of each haplotype state
    @return: None
    """
    # create multinomial probabilities using counts and fitnesses
    nmutations = len(mutations)
    probs = np.zeros(nmutations)
    for i, m in enumerate(mutations):
        probs[i] = m.count * fitnesses[m.state]
    probs /= probs.sum()
    # do multinomial sampling corresponding to new counts for the mutations
    counts = np.random.multinomial(N_hap, probs)
    # update the mutation counts
    for i in range(nmutations):
        mutations[i].count = counts[i]
    return

def _deep_prune(dead_line):
    """
    @param dead_line: a dead mutational line
    """
    while True:
        parent = dead_line.parent
        if parent is None:
            raise ValueError('pruning the original background mutational line')
        if id(dead_line) in parent.child_id_map:
            del parent.child_id_map[id(dead_line)]
            if parent.is_dead():
                dead_line = parent
            else:
                return
        else:
            return

def prune(mutations_in):
    """
    Prune mutational lines which are dead by descent.
    @param mutations: a sequence of mutation objects to be updated
    @return: sequence of extant mutational lines
    """
    mutations_out = []
    for m in mutations_in:
        if m.is_dead():
            _deep_prune(m)
        else:
            mutations_out.append(m)
    return mutations_out

def evolve(N_hap, mu, generation, mutations, fitnesses, fitness_history):
    """
    Evolve one generation.
    @param N_hap: total haplotype population size
    @param mu: expected total mutations per generation
    @param generation: the current generation
    @param mutations: a sequence of mutation objects
    @param fitnesses: fitness of each haplotype state
    @param fitness_history: tracks the average fitness at each generation
    @return: a new list of extant mutational lines by descent
    """
    mutate(mu, generation, mutations)
    reselect(N_hap, mutations, fitnesses)
    mutations = prune(mutations)
    # update the fitness history
    total_fitness = sum(m.count * fitnesses[m.state] for m in mutations)
    fitness_history.append(total_fitness / N_hap)
    # return the sequence of new mutational lines
    return mutations

class TestWrightFisherForwardKLineage(unittest.TestCase):
    def test_forward(self):
        pass

if __name__ == '__main__':
    unittest.main()

