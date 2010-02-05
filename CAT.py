"""
This module is written to test some hypothesis about the CAT model.
"""

from StringIO import StringIO
import math

import numpy as np

import Util

g_filename = '/home/argriffi/packages/phylobayes2.3f/sources/C20.ProfilesCow72'
g_aa_letters_a = 'ARNDCQEGHILKMFPSTWYV'
g_aa_letters_b = 'ACDEFGHIKLMNPQRSTVWY'

def almost_equals(a, b, eps=1e-8):
    return abs(b-a) < eps

def assert_profile(profile, eps=1e-8):
    profile_sum = sum(profile)
    assert almost_equals(profile_sum, 1, eps), profile_sum

def profile_power(profile, exponent):
    """
    @param profile: a profile as a numpy array
    @param exponent: the power to which to raise each element
    """
    raw_exp_profile = profile**exponent
    exp_profile = raw_exp_profile / sum(raw_exp_profile)
    assert_profile(exp_profile)
    return exp_profile

def main():
    fin = open(g_filename)
    lines = list(Util.stripped_lines(fin.readlines()))
    float_rows = [[float(x) for x in line.split()] for line in lines[1:]]
    mixture_profile = np.array([row[0] for row in float_rows])
    aa_profiles = np.array([row[1:] for row in float_rows])
    # assert that we have the right number of elements
    assert mixture_profile.shape == (20,)
    assert aa_profiles.shape == (20, 20)
    # assert that profiles sum to one
    assert_profile(mixture_profile)
    for aa_profile in aa_profiles:
        assert_profile(aa_profile)
    for exponent in [0.5, 1.0, 2.0, 10.0]:
        print 'exponent', exponent
        print
        # get the total profile
        total_profile = np.zeros(20)
        for weight, profile in zip(mixture_profile, aa_profiles):
            total_profile += weight * profile
        # get the exponentiated total profile
        exp_total_profile = profile_power(total_profile, exponent)
        # get the weighted sum of exponentiated component profiles
        exp_component_profile = np.zeros(20)
        for weight, profile in zip(mixture_profile, aa_profiles):
            exp_component_profile += weight * profile_power(profile, exponent)
        # assert that we have two valid profiles
        assert_profile(total_profile)
        assert_profile(exp_total_profile)
        assert_profile(exp_component_profile)
        # compare the profiles elementwise by log odds
        elementwise_log_odds = [math.log(b/a) for a, b in zip(total_profile, exp_component_profile)]
        for aa_letter, log_odds in zip(g_aa_letters_b, elementwise_log_odds):
            print aa_letter, log_odds
        print
        for aa_letter, log_odds in zip(g_aa_letters_b, elementwise_log_odds):
            if aa_letter in 'IVYWREL':
                print aa_letter, log_odds
        print

if __name__ == '__main__':
    main()
