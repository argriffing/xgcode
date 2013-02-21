"""
Check properties of a rate matrix augmented with hidden binary variables.

The binary variables allow or disallow transitions to corresponding states
of non-hidden variables.
For example, this could be used when some states are invisibly
disallowed at various times for whatever reason,
and the set of disallowed states changes
according to a continuous time process,
and only one state changes allowance status at any instant.
If the original rate matrix is Q with N states,
then the rate matrix augmented with these hidden binary variables
will be called Q* and it will have N * 2^(N-1) states.
The coefficient N corresponds to the N states of the visible variable,
and the term 2^(N-1) corresponds to the possible binary allowed/disallowed
states of the binary variables corresponding to the other N-1 states
of the visible variable.
The exponent is N-1 instead of N because the current state of the
visible variable cannot be disallowed.
So at most N-1 states of the visible variable can be disallowed
at any given time.
"""

import itertools

import numpy as np

import MatrixUtil


def get_expanded_states(n):
    """
    @param n: number of states
    @return: a permutation and a binary vector for each state
    """
    perms = list(itertools.permutations(range(n)))
    masks = [(1,) + mask for mask in itertools.product((0,1), repeat=n-1)]
    return tuple(itertools.product(perms, masks))

def get_toggle_on_adjacency(expo_states):
    """
    Precompute this design matrix.
    @param expo_states: expanded states with permutations and masks
    @return: a binary adjacency matrix
    """
    nstates = len(expo_states)
    adjacency = np.zeros((nstates, nstates), dtype=int)
    for i in range(nstates):
        for j in range(nstates):
            perm_i, mask_i = expo_states[i]
            perm_j, mask_j = expo_states[j]
            ntog_10 = sum(1 for a, b in zip(mask_i, mask_j) if a==1 and b==0)
            ntog_01 = sum(1 for a, b in zip(mask_i, mask_j) if a==0 and b==1)
            if perm_i == perm_j and ntog_01 == 1 and ntog_10 == 0:
                adjacency[i, j] = 1
    return adjacency

def get_toggle_off_adjacency(expo_states):
    """
    Precompute this design matrix.
    @param expo_states: expanded states with permutations and masks
    @return: a binary adjacency matrix
    """
    nstates = len(expo_states)
    adjacency = np.zeros((nstates, nstates), dtype=int)
    for i in range(nstates):
        for j in range(nstates):
            perm_i, mask_i = expo_states[i]
            perm_j, mask_j = expo_states[j]
            ntog_10 = sum(1 for a, b in zip(mask_i, mask_j) if a==1 and b==0)
            ntog_01 = sum(1 for a, b in zip(mask_i, mask_j) if a==0 and b==1)
            if perm_i == perm_j and ntog_01 == 0 and ntog_10 == 1:
                adjacency[i, j] = 1
    return adjacency

def get_observable_adjacency(expo_states):
    """
    Precompute this design matrix.
    Entries of the matrix will be 1 when the visible state changes.
    @param expo_states: expanded states with permutations and masks
    @return: a binary adjacency matrix
    """
    nstates = len(expo_states)
    adjacency = np.zeros((nstates, nstates), dtype=int)
    for i in range(nstates):
        for j in range(nstates):
            perm_i, mask_i = expo_states[i]
            perm_j, mask_j = expo_states[j]
            if mask_i != mask_j:
                continue
            if perm_i[0] != perm_j[0] and hamdist(perm_i, perm_j) == 2:
                nvis = len(perm_i)
                kflip = None
                for k in range(nvis):
                    if (perm_i[0], perm_i[k]) == (perm_j[k], perm_j[0]):
                        kflip = k
                if kflip is None:
                    raise Exception
                if mask_i[kflip] == 1 and mask_j[kflip] == 1:
                    adjacency[i, j] = 1
    return adjacency

def hamdist(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)

def blinkize(pre_Q, blink_birth, blink_death):
    """
    @param pre_Q: a pre-rate matrix
    @param blink_birth: a rate
    @param blink_death: a rate
    @return: a new 'blink-ized' pre-rate matrix
    """
    n = pre_Q.shape[0]
    
    # pre-compute some properties of the new rate matrix
    expo_states = get_expanded_states(n)
    adj_toggle_on = get_toggle_on_adjacency(expo_states)
    adj_toggle_off = get_toggle_off_adjacency(expo_states)
    adj_observable = get_observable_adjacency(expo_states)

    # check some invariants
    if np.any(adj_toggle_on * adj_toggle_off):
        raise Exception
    adj_toggle_either = adj_toggle_on + adj_toggle_off
    MatrixUtil.assert_symmetric(adj_toggle_either)
    MatrixUtil.assert_symmetric(adj_observable)
    if np.any(adj_toggle_either * adj_observable):
        raise Exception

    # initialize the new rate matrix
    nblinks = len(expo_states)
    pre_blink = np.zeros((nblinks, nblinks), dtype=float)
    
    # add the rates that control the toggling between the hidden states
    pre_blink += adj_toggle_off * blink_death
    pre_blink += adj_toggle_on * blink_birth

    # add the rates that control the observable state change
    for i in range(nblinks):
        for j in range(nblinks):
            perm_i, mask_i = expo_states[i]
            perm_j, mask_j = expo_states[j]
            if adj_observable[i, j]:
                a = perm_i[0]
                b = perm_j[0]
                pre_blink[i, j] = pre_Q[a, b]

    # return the newly constructed pre-rate matrix
    return pre_blink

def pre_Q_to_Q(pre_Q):
    return pre_Q - np.diag(np.sum(pre_Q, axis=1))

def sample_reversible_pre_Q(n):
    M = np.exp(np.random.randn(n, n))
    S = M + M.T
    weights = np.exp(np.random.randn(n))
    v = weights / np.sum(weights)
    pre_Q = np.dot(S, np.diag(v))
    return pre_Q

