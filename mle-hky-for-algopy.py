"""
This is an ML estimate of HKY85 molecular evolutionary parameters.
"""

import numpy as np
from scipy import optimize, linalg

# AGCT subsititution counts between human and chimp mitochondrial coding dna.
g_data = [
        [2954, 141, 17, 16],
        [165, 1110, 5, 2],
        [18, 4, 3163, 374],
        [15, 2, 310, 2411],
        ]

def transform_params(Y):
    X = np.exp(Y)
    tsrate, tvrate = X[0], X[1]
    v_unnormalized = np.array([X[2], X[3], X[4], 1.0])
    v = v_unnormalized / np.sum(v_unnormalized)
    return tsrate, tvrate, v

def minimize_me(Y):
    a, b, v = transform_params(Y)
    Q = np.array([
        [0, a, b, b],
        [a, 0, b, b],
        [b, b, 0, a],
        [b, b, a, 0],
        ])
    Q = np.dot(Q, np.diag(v))
    Q -= np.diag(np.sum(Q, axis=1))
    S = np.log(np.dot(np.diag(v), linalg.expm(Q)))
    return -np.sum(S * g_data)

def main():
    results = optimize.fmin(
            minimize_me, np.zeros(5),
            maxiter=10000, maxfun=10000, full_output=True)
    tsrate, tvrate, v = transform_params(results[0])
    print 'results output from fmin:', results
    print 'estimated transition rate parameter:', tsrate
    print 'estimated transversion rate parameter:', tvrate
    print 'estimated stationary distribution:', v

if __name__ == '__main__':
    main()

