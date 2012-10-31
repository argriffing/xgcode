"""
Construct a nucleotide rate matrix with recessivity-inspired parameters.

This model of nucleotide evolution is traditional in the sense that it is
a time-reversible 4-state continuous-time Markov process,
and its equilibrium state distribution corresponds
directly to a subset of the model parameters.
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
drift, natural selection, dominance/recessivity, migration,
population size, generation time, etc.
but for now I am looking at only a 4-state model of nucleotide evolution
with only a couple of parameters that do not correspond directly to
the overall substitution rate or to the equilibrium nucleotide distribution.
All seven options listed below are inconveniently transformed in a way that
allows any positive or negative number to be used for any option.
The d parameter is better known as 2*h - 1, and is defined with respect
to a new mutant that is more fit than the previously fixed allele.
Note that if h' = 1-h then d' = -d.
The k parameter controls the sharpness of the recessivity/dominance flip
as a selection parameter changes sign;
when k is a large positive number this change is very abrupt,
whereas a negative value of k with large absolute value means that
selection is nearly genic for a wider range of selection values.
When the log expected rate is zero,
then the rate matrix will be scaled such that the process has
one expected nucleotide substitution per time unit.
When either d=0 or k<<0 this model should reduce to the
rate matrix constructed using the neutral Kimura-inspired parameterization
under a Jukes-Cantor mutation model.
"""

from StringIO import StringIO
import math

import numpy
import scipy
import scipy.special

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
            Form.Float('d', 'dominance d of fitter allele', 0.4),
            Form.Float('k', 'Kacser and Burn effect parameter k', 2.0),
            Form.Float('log_rate', 'log expected rate', 0),
            ]

def get_form_out():
    return FormOut.Report()

def get_fixation_unconstrained(S, d, k):
    D = d * numpy.tanh(numpy.exp(k)*S)
    H = numpy.zeros_like(S)
    for i in range(H.shape[0]):
        for j in range(H.shape[1]):
            H[i, j] = 1. / kimrecessive.denom_piecewise(
                    0.5*S[i, j], D[i, j])
    return H

def create_pre_rate_matrix(fs):
    """
    The returned rate matrix is not diagonally adjusted or rate-scaled.
    """

    # Construct the stationary distribution
    log_v = numpy.array([
        fs.A_log_weight,
        fs.C_log_weight,
        fs.G_log_weight,
        fs.T_log_weight,
        ])
    v = numpy.exp(log_v)
    v /= numpy.sum(v)
    
    # Construct F and S in the Yang-Nielsen 2008 notation.
    F = log_v
    e = numpy.ones_like(F)
    S = numpy.outer(e, F) - numpy.outer(F, e)

    # Compute the recessivity component of the rate matrix.
    H = get_fixation_unconstrained(S, fs.d, fs.k)

    # Construct and return the pre rate matrix and stationary distribution.
    pre_Q = H
    return H, v


def get_response_content(fs):
    numpy.set_printoptions(linewidth=200)
    out = StringIO()
    #
    pre_Q, v = create_pre_rate_matrix(fs)
    Q = pre_Q - numpy.diag(numpy.sum(pre_Q, axis=1))
    observed_rate = -numpy.dot(v, numpy.diag(Q))
    requested_rate = numpy.exp(fs.log_rate)
    Q *= (requested_rate / observed_rate)
    #
    print >> out, 'stationary distribution:'
    print >> out, v
    print >> out
    print >> out, 'rate matrix:'
    print >> out, Q
    print >> out
    print >> out, 'pi_i Q_ij (to check reversibility):'
    print >> out, numpy.dot(numpy.diag(v), Q)
    print >> out
    """
    print >> out, '<html><head></head><body>'
    print >> out, '<pre>'
    print >> out, 'stationary distribution:'
    print >> out, v
    print >> out, '</pre>'
    print >> out
    print >> out, '<table>'
    print >> out, '</table>'
    print >> out, '</body></html>'
    """
    #
    return out.getvalue()

