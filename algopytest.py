"""
The parameter names
a, c, e, l, b, d, k, s
are an attempt to translate the greek letters
in the Schadt et al. paper.
doi:10.1101/gr.8.3.222
"""

import numpy as np

import sympy


# define some indices
g_A, g_G, g_C, g_T = 0, 1, 2, 3

# These are human->chimp AGCT mitochondrial sequence nucleotide differences.
# We should really use sequence data from the tips of a rooted tree
# and use Felsenstein's pruning algorithm to integrate over all
# possible combinations of ancestral states,
# but that would be too complicated for our purposes here.
g_data = [
        [2954, 141, 17, 16],
        [165, 1110, 5, 2],
        [18, 4, 3163, 374],
        [15, 2, 310, 2411],
        ]

def get_stationary_distribution(X):
    a, c, e, l, b, d, k, s = X.tolist()
    denom_AG = (a + c + e + l) * (c + d + k + l)
    denom_CT = (b + d + k + s) * (c + d + k + l)
    pi_A = ( e * (d + k) + d * (c + l) ) / denom_AG
    pi_G = ( a * (d + k) + k * (c + l) ) / denom_AG
    pi_C = ( c * (d + k) + s * (c + l) ) / denom_CT
    pi_T = ( l * (d + k) + b * (c + l) ) / denom_CT
    return np.array([pi_A, pi_G, pi_C, pi_T])

def get_P(X):
    """
    Manually type the 16 entries of the P matrix.
    """
    t = 1
    P = np.zeros((4, 4))
    a, c, e, l, b, d, k, s = X.tolist()
    pi_A, pi_G, pi_C, pi_T = get_stationary_distribution(X).tolist()
    c1 = d + k + c + l
    c2 = c + l
    c3 = a + c + e + l
    c4 = k - a
    c5 = d + k
    c6 = d + k + s + b
    c7 = c - s
    c8 = d - e
    c9 = l - b
    P[g_A, g_A] = sum(
            pi_A,
            -( (c2*c8) / (c1*(c3-c1)) ) * exp(-c1*t),
            ( (a*(c3 - c1) - c2*c4) / (c3*(c3-c1)) ) * exp(-c3*t),
            )
    P[g_A, g_G] = sum(
            pi_G,
            -( (c2*c4) / (c1*(c3-c1)) ) * exp(-c1*t),
            ( (c2*c4 - a*(c3 - c1)) / (c3*(c3-c1)) ) * exp(-c3*t),
            )
    return P

def check_sympy_derivative():
    t = 1
    a, c, e, l, b, d, k, s = sympy.symbols('a c e l b d k s')
    denom_AG = (a + c + e + l) * (c + d + k + l)
    denom_CT = (b + d + k + s) * (c + d + k + l)
    pi_A = ( e * (d + k) + d * (c + l) ) / denom_AG
    pi_G = ( a * (d + k) + k * (c + l) ) / denom_AG
    pi_C = ( c * (d + k) + s * (c + l) ) / denom_CT
    pi_T = ( l * (d + k) + b * (c + l) ) / denom_CT
    c1 = d + k + c + l
    c2 = c + l
    c3 = a + c + e + l
    c4 = k - a
    c5 = d + k
    c6 = d + k + s + b
    c7 = c - s
    c8 = d - e
    c9 = l - b
    p_AG = pi_G - ( (c2*c4) / (c1*(c3-c1)) ) * sympy.exp(-c1*t) + (
            ( (c2*c4 - a*(c3 - c1)) / (c3*(c3-c1)) ) * sympy.exp(-c3*t) )
    my_d = sympy.diff(p_AG, a)
    #a, c, e, l, b, d, k, s = 1, 1, 1, 1, 1, 1, 1, 1
    #print my_d
    #print sympy.N(my_d)
    #print sympy.evalf(a+c, a=1, c=2)
    #print (a+c).evalf(a=1, c=2)
    #print (a+c).evalf()
    #print (a+c).subs((a, c), (1, 2))
    #print (a+c).subs(a, 1).subs(c, 2)

def main():
    check_sympy_derivative()

if __name__ == '__main__':
    main()

