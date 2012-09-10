"""
Check W-F compensatory pathway preference as a function of population size.

Do not use forward simulation.
Presumably the exact solutions have been checked against forward simulation
in other web scripts.
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
import wfbckcompens
import combobreaker

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('theta', 'mutation rate 4*N*mu', '1.0',
                low_inclusive=0.001),
            Form.Float('Nr', 'recombination rate N*r', '5.0',
                low_inclusive=0),
            Form.Float('Ns', 'selection N*s', '1.0',
                low_inclusive=0),
            ]

def get_form_out():
    return FormOut.Report()

def get_p2(N_diploid, theta, Nr, Ns):
    """
    Compute the probability of compensatory substitution pathway type 2.
    @param N_diploid: diploid population size
    @param theta: a mutation rate
    @param Nr: a recombination rate
    @param Ns: a selection value
    @return: p_t2
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
    X = linalg.solve(np.eye(nstates - k) - P[k:, k:], P[k:, :k])
    H = P[:k, :k] + np.dot(P[:k, k:], X)
    p_direct = H[0, 3] / (1 - H[0,0])
    # The following line is Equation (1) of the Nasrallah manuscript.
    p_t2 = (2*p_direct) / (1 + p_direct)
    p_t1 = 1 - p_t2
    #
    return p_t2

def get_response_content(fs):
    out = StringIO()
    print >> out, 'theta:', fs.theta
    print >> out, 'Nr:', fs.Nr
    print >> out, 'Ns:', fs.Ns
    print >> out
    for N_diploid in (3, 4, 5, 6, 7, 8, 9, 10, 11, 12):
        N_hap = N_diploid * 2
        # TODO check s vs N_diploid so that selection is not too big
        # as to make fitnesses negative or zero
        #
        p2 = get_p2(N_diploid, fs.theta, fs.Nr, fs.Ns)
        #
        print >> out, '2N:', N_hap
        print >> out, 'Type 2 probability:', p2
        print >> out
    return out.getvalue()


def main():
    """
    Check some stuff that takes longer to compute.
    """
    pass

if __name__ == '__main__':
    main()

