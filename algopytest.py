"""
The parameter names
a, c, e, l, b, d, k, s
are an attempt to translate the greek letters
in the Schadt et al. paper.
doi:10.1101/gr.8.3.222
"""

import numpy as np

import sympy
import theano
import theano.tensor as T
from theano import pp


# define some indices
g_A, g_G, g_C, g_T = 0, 1, 2, 3

# These are human->chimp AGCT mitochondrial sequence nucleotide differences.
# We should really use sequence data from the tips of a rooted tree
# and use Felsenstein's pruning algorithm to integrate over all
# possible combinations of ancestral states,
# but that would be too complicated for our purposes here.
g_data = [
        2954, 141, 17, 16,
        165, 1110, 5, 2,
        18, 4, 3163, 374,
        15, 2, 310, 2411,
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

def check_theano_gradient():
    a = T.dscalar('a')
    c = T.dscalar('c')
    e = T.dscalar('e')
    l = T.dscalar('l')
    b = T.dscalar('b')
    d = T.dscalar('d')
    k = T.dscalar('k')
    s = T.dscalar('s')
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
    #
    # numerators of a part of the transition matrix
    p_AA_n1 = c2*c8
    p_AG_n1 = c2*c4
    p_AC_n1 = c2*c7
    p_AT_n1 = c2*c9
    p_GA_n1 = c2*c8
    p_GG_n1 = c2*c4
    p_GC_n1 = c2*c7
    p_GT_n1 = c2*c9
    p_CA_n1 = c5*c8
    p_CG_n1 = c5*c4
    p_CC_n1 = c5*c7
    p_CT_n1 = c5*c9
    p_TA_n1 = c5*c8
    p_TG_n1 = c5*c4
    p_TC_n1 = c5*c7
    p_TT_n1 = c5*c9
    #
    # denominators of a part of the transition matrix
    p_AA_d1 = c1 * (c3 - c1)
    p_AG_d1 = c1 * (c3 - c1)
    p_AC_d1 = c1 * (c6 - c1)
    p_AT_d1 = c1 * (c6 - c1)
    p_GA_d1 = c1 * (c3 - c1)
    p_GG_d1 = c1 * (c3 - c1)
    p_GC_d1 = c1 * (c6 - c1)
    p_GT_d1 = c1 * (c6 - c1)
    p_CA_d1 = c1 * (c3 - c1)
    p_CG_d1 = c1 * (c3 - c1)
    p_CC_d1 = c1 * (c6 - c1)
    p_CT_d1 = c1 * (c6 - c1)
    p_TA_d1 = c1 * (c3 - c1)
    p_TG_d1 = c1 * (c3 - c1)
    p_TC_d1 = c1 * (c6 - c1)
    p_TT_d1 = c1 * (c6 - c1)
    #
    # numerators of a different part of the transition matrix
    p_AA_n2 = a*(c3 - c1) - c2*c4
    p_AG_n2 = c2*c4 - a*(c3 - c1)
    p_AC_n2 = l*s - c*b
    p_AT_n2 = b*c - l*s
    p_GA_n2 = c2*c8 - e*(c3 - c1)
    p_GG_n2 = e*(c3 - c1) - c2*c8
    p_GC_n2 = l*s - c*b
    p_GT_n2 = b*c - l*s
    p_CA_n2 = e*k - d*a
    p_CG_n2 = d*a - e*k
    p_CC_n2 = b*(c6 - c1) - c5*c9
    p_CT_n2 = c5*c9 - b*(c6 - c1) # typo in paper c7 should be c6 ?
    p_TA_n2 = e*k - d*a
    p_TG_n2 = d*a - e*k
    p_TC_n2 = c5*c7 - s*(c6 - c1)
    p_TT_n2 = s*(c6 - c1) - c5*c7
    #
    # denominators of a different part of the transition matrix
    p_AA_d2 = c3 * (c3 - c1)
    p_AG_d2 = c3 * (c3 - c1)
    p_AC_d2 = c6 * (c6 - c1)
    p_AT_d2 = c6 * (c6 - c1)
    p_GA_d2 = c3 * (c3 - c1)
    p_GG_d2 = c3 * (c3 - c1)
    p_GC_d2 = c6 * (c6 - c1)
    p_GT_d2 = c6 * (c6 - c1)
    p_CA_d2 = c3 * (c3 - c1)
    p_CG_d2 = c3 * (c3 - c1)
    p_CC_d2 = c6 * (c6 - c1)
    p_CT_d2 = c6 * (c6 - c1)
    p_TA_d2 = c3 * (c3 - c1)
    p_TG_d2 = c3 * (c3 - c1)
    p_TC_d2 = c6 * (c6 - c1)
    p_TT_d2 = c6 * (c6 - c1)
    #
    # probabilities
    p_AA = pi_A - (p_AA_n1/p_AA_d1) * T.exp(-c1) + (p_AA_n2/p_AA_d2)*T.exp(-c3)
    p_AG = pi_G - (p_AG_n1/p_AG_d1) * T.exp(-c1) + (p_AG_n2/p_AG_d2)*T.exp(-c3)
    p_AC = pi_C + (p_AC_n1/p_AC_d1) * T.exp(-c1) + (p_AC_n2/p_AC_d2)*T.exp(-c6)
    p_AT = pi_T + (p_AT_n1/p_AT_d1) * T.exp(-c1) + (p_AT_n2/p_AT_d2)*T.exp(-c6)
    p_GA = pi_A - (p_GA_n1/p_GA_d1) * T.exp(-c1) + (p_GA_n2/p_GA_d2)*T.exp(-c3)
    p_GG = pi_G - (p_GG_n1/p_GG_d1) * T.exp(-c1) + (p_GG_n2/p_GG_d2)*T.exp(-c3)
    p_GC = pi_C + (p_GC_n1/p_GC_d1) * T.exp(-c1) + (p_GC_n2/p_GC_d2)*T.exp(-c6)
    p_GT = pi_T + (p_GT_n1/p_GT_d1) * T.exp(-c1) + (p_GT_n2/p_GT_d2)*T.exp(-c6)
    p_CA = pi_A - (p_CA_n1/p_CA_d1) * T.exp(-c1) + (p_CA_n2/p_CA_d2)*T.exp(-c3)
    p_CG = pi_G - (p_CG_n1/p_CG_d1) * T.exp(-c1) + (p_CG_n2/p_CG_d2)*T.exp(-c3)
    p_CC = pi_C + (p_CC_n1/p_CC_d1) * T.exp(-c1) + (p_CC_n2/p_CC_d2)*T.exp(-c6)
    p_CT = pi_T + (p_CT_n1/p_CT_d1) * T.exp(-c1) + (p_CT_n2/p_CT_d2)*T.exp(-c6)
    p_TA = pi_A - (p_TA_n1/p_TA_d1) * T.exp(-c1) + (p_TA_n2/p_TA_d2)*T.exp(-c3)
    p_TG = pi_G - (p_TG_n1/p_TG_d1) * T.exp(-c1) + (p_TG_n2/p_TG_d2)*T.exp(-c3)
    p_TC = pi_C + (p_TC_n1/p_TC_d1) * T.exp(-c1) + (p_TC_n2/p_TC_d2)*T.exp(-c6)
    p_TT = pi_T + (p_TT_n1/p_TT_d1) * T.exp(-c1) + (p_TT_n2/p_TT_d2)*T.exp(-c6)
    #
    #
    print pp(p_AG)
    mydict = {
        a : 1.0,
        c : 1.0,
        e : 1.0,
        l : 1.0,
        b : 1.0,
        d : 1.0,
        k : 1.0,
        s : 1.0,
        }
    print T.grad(p_AG).eval({
        a : 1.0,
        c : 1.1,
        e : 1.2,
        l : 1.3,
        #b : 1.0,
        d : 1.4,
        k : 1.5,
        #s : 1.0,
        })
    #print p_AG.eval(mydict)
    #
    # check the stationary distribution when all rates are the same
    print pi_A.eval({
        a : 1.0,
        c : 1.0,
        e : 1.0,
        l : 1.0,
        #b : 1.0,
        d : 1.0,
        k : 1.0,
        #s : 1.0,
        })
    print pi_G.eval({
        a : 1.0,
        c : 1.0,
        e : 1.0,
        l : 1.0,
        #b : 1.0,
        d : 1.0,
        k : 1.0,
        #s : 1.0,
        })
    print pi_C.eval({
        #a : 1.0,
        c : 1.0,
        #e : 1.0,
        l : 1.0,
        b : 1.0,
        d : 1.0,
        k : 1.0,
        s : 1.0,
        })
    print pi_T.eval({
        #a : 1.0,
        c : 1.0,
        #e : 1.0,
        l : 1.0,
        b : 1.0,
        d : 1.0,
        k : 1.0,
        s : 1.0,
        })
    #
    # check that denominators are not zero
    print p_AG_d1.eval({
        a : 1.0,
        c : 1.0,
        e : 1.0,
        l : 1.0,
        #b : 1.0,
        d : 1.0,
        k : 1.0,
        #s : 1.0,
        })
    print p_AC_d1.eval({
        #a : 1.0,
        c : 1.0,
        #e : 1.0,
        l : 1.0,
        b : 1.0,
        d : 1.0,
        k : 1.0,
        s : 1.0,
        })
    print p_AG_d2.eval({
        a : 1.0,
        c : 1.0,
        e : 1.0,
        l : 1.0,
        #b : 1.0,
        d : 1.0,
        k : 1.0,
        #s : 1.0,
        })
    print p_AC_d2.eval({
        #a : 1.0,
        c : 1.0,
        #e : 1.0,
        l : 1.0,
        b : 1.0,
        d : 1.0,
        k : 1.0,
        s : 1.0,
        })
    """
    my_ll = T.dot(
            [
                foo,
                bar,
                ],
            g_data,
            )
    """

    


def main():
    #check_sympy_derivative()
    check_theano_gradient()

if __name__ == '__main__':
    main()

