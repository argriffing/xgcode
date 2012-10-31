"""
Construct a nucleotide rate matrix with recessivity-inspired parameters. [UNFINISHED]

This model of nucleotide evolution is traditional in the sense that it is
a time-reversible 4-state continuous-time Markov process,
and its equilibrium state distribution corresponds
directly to a subset of model parameters.
It is also traditional in the sense that there is a parameter that
can be adjusted to simultaneously rescale all of the rates by the same amount.
But this model also has two parameters which are inspired by
the population genetic theoretic fixation probability of new mutants
in the presence of selection with some dominance or recessivity
of the more fit allele with respect to the less fit allele.
This diffusion theory approximation of Wright-Fisher population genetic models
is explained in a 1957 paper by Kimura.
Perhaps eventually these kinds of parameters may make their way into
more complicated phylogenetic models of evolution that try to model
the interaction of many population-genetic concepts like mutation, fixation,
drift, selection, dominance/recessivity, migration, population size, etc.
but for now I am looking at only a 4-state model of nucleotide evolution
with only a couple of parameters that do not correspond directly to
the overall substitution rate or to the equilibrium nucleotide distribution.
All seven listed parameters are inconveniently transformed in a way that
allows any positive or negative number to be used for any parameter.
The d parameter is better known as 2*h - 1, and is defined with respect
to a new mutant that is more fit than the previously fixed allele.
The k parameter controls the sharpness of the recessivity/dominance flip
as a selection parameter changes sign;
when k is a large positive number this change is very abrupt,
whereas the change is nearly linear when k is a negative number with large
absolute value.
"""

from StringIO import StringIO
import math
import cmath
import random

import numpy
import scipy
import scipy.special
import algopy

import Form
import FormOut
import kimrecessive

def get_form():
    """
    @return: the body of a form
    """
    return [
            Form.Float('A_log_weight', 'A log weight', -1.1),
            Form.Float('C_log_weight', 'C log weight', 1.2),
            Form.Float('G_log_weight', 'G log weight', 0.0),
            Form.Float('T_log_weight', 'T log weight', 0.3),
            Form.Float('d', 'dominance of fitter allele d', 0.4),
            Form.Float('k', 'Kacser and Burn effect parameter k', 2.0),
            Form.Float('log_rate', 'log expected rate', 0),
            ]

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    #numpy.set_printoptions(linewidth=200)
    out = StringIO()
    #
    print >> out, 'generating random selection and dominance...'
    print >> out
    for i in range(10):
        #c = random.random()
        #d = random.random()
        c = random.gauss(0, 1)
        d = random.gauss(0, 0.01)
        print >> out, 'c:', c
        print >> out, 'd:', d
        print >> out
        print >> out, 'method a:'
        print >> out, kimura_1957_54_denominator_analytic(c, d)
        print >> out
        print >> out, 'method b:'
        print >> out, kimura_1957_54_denominator_analytic_b(c, d)
        print >> out
        print >> out, 'asymptotic:'
        print >> out, kimura_1957_54_denominator_analytic_asymptotic(c, d)
        print >> out
        print >> out, 'asymptotic, genic:'
        print >> out, kimura_1957_54_denominator_analytic_asymptotic(c, 0)
        print >> out
        print >> out, 'explicitly genic:'
        print >> out, kimura_genic(c)
        print >> out
        print >> out
    print >> out, 'asymptotic(c=0.5, d=0):'
    print >> out, kimura_1957_54_denominator_analytic_asymptotic(0.5, 0)
    print >> out
    #print >> out, 'general_purpose(c=1, d=0):'
    #print >> out, kimura_1957_54_denominator_analytic_b(1, 0)
    #print >> out, '(division by zero error)'
    #print >> out
    print >> out, 'general_purpose(c=0, d=0.5):'
    print >> out, kimura_1957_54_denominator_analytic_b(0, 0.5)
    print >> out
    #
    return out.getvalue()

