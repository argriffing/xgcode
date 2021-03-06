"""
Analyze multiple sequential nucleotide sites of a pileup.

Parameters x, y, and z are branch lengths which define the two three-taxon
trees whose Jukes-Cantor nucleotide distribution at the tips
define the {RR, RA, AA, AB} zygosity distributions for the
recent vs. ancient mrca states.
The low, medium, and high parameters are expectations of three geometrically
distributed mixture components of a garbage state.
The seqerror parameter is the probability of sequencing randomization;
this the probability that the sequencing machine spits out a random
nucleotide instead of the correct nucleotide.
The nomcoverage parameter defines the nominal coverage of the pileup.
The kmulticoverages parameter defines the number of
nominal coverage multiples which might result from duplications.
"""

from StringIO import StringIO
import itertools
import math
import argparse

import scipy.misc

from SnippetUtil import HandlingError
import Form
import FormOut
import DGRP
import ReadCoverageRef
import Util
import ExternalHMM
import lineario
import TransitionMatrix

g_default_params = [
        ('x', '0.1'),
        ('y', '0.01'),
        ('z', '0.0001'),
        ('low', '2'),
        ('med', '20'),
        ('high', '1000'),
        ('seqerr', '.01'),
        ('nomcoverage', '20'),
        ('kmulticoverages', '4')]

g_default_data = [
        (20, 0, 0, 0),
        (22, 2, 2, 1),
        (20, 1, 0, 5),
        (6, 14, 0, 0),
        (21, 0, 2, 1),
        (28, 1, 2, 1),
        (21, 2, 0, 1),
        (20, 0, 0, 0),
        (12, 12, 0, 1),
        (20, 1, 0, 5),
        (6, 14, 0, 0),
        (9, 0, 11, 1),
        (28, 1, 2, 1),
        (21, 2, 0, 1),
        (20, 0, 0, 0),
        (12, 12, 5, 1),
        (20, 30, 0, 5),
        (6, 14, 0, 0),
        (9, 0, 0, 1),
        (4, 1, 2, 1),
        (21, 2, 0, 1)]

def get_form():
    """
    @return: the body of a form
    """
    data_lines = ['\t'.join(str(x) for x in row) for row in g_default_data]
    form_objects = [
            Form.Integer('region_size', 'expected region size',
                10, low=1, high=1000000),
            Form.Float('misalignment_effect', 'misalignment effect',
                '0.5', low_inclusive=0),
            Form.MultiLine('param_field', 'parameters',
                '\n'.join('\t'.join(p) for p in g_default_params)),
            Form.MultiLine('data_field', 'data',
                '\n'.join(data_lines))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    lines = Util.get_stripped_lines(StringIO(fs.param_field))
    model = DGRP.Model()
    model.from_lines(lines)
    # see how the states interact with the observations
    states = (
            model.get_recent_state(),
            model.get_ancient_state(),
            model.get_misaligned_state(fs.misalignment_effect),
            model.get_garbage_state())
    # define the transition object
    nstates = len(states)
    prandom = min(1.0, (nstates / (nstates - 1.0)) / fs.region_size)
    T = TransitionMatrix.UniformTransitionObject(prandom, nstates)
    # use StringIO objects for storage
    hmm = ExternalHMM.ExternalModel(T, states, (None, None, None))
    converter = lineario.IntTupleConverter()
    o_stream = lineario.SequentialStringIO(converter, fs.data_field)
    hmm.init_dp(o_stream)
    o_stream.open_read()
    for p, obs in itertools.izip(hmm.posterior(), o_stream.read_forward()):
        p_recent, p_ancient, p_misaligned, p_garbage = p
        # get the prior probability of polymorphism conditional on state
        p_recent_AA = states[0].get_posterior_distribution(obs)[2]
        p_ancient_AA = states[1].get_posterior_distribution(obs)[2]
        # compute the posterior probability of a polymorphism
        posterior_polymorphism = 0
        posterior_polymorphism += p_recent * p_recent_AA
        posterior_polymorphism += p_ancient * p_ancient_AA
        # Given that a polymorphism occurred,
        # get the probability distribution over the
        # three non-reference nucleotides.
        r = model.seqerr
        log_Pr = math.log(r/4.0)
        log_PA = math.log(1 - 3*r/4.0)
        logs = [
                obs[1]*log_PA + obs[2]*log_Pr + obs[3]*log_Pr,
                obs[1]*log_Pr + obs[2]*log_PA + obs[3]*log_Pr,
                obs[1]*log_Pr + obs[2]*log_Pr + obs[3]*log_PA]
        condmaxpost = math.exp(max(logs) - scipy.misc.logsumexp(logs))
        # get the posterior probability distribution
        maxpost = posterior_polymorphism * condmaxpost
        # show the inference for this position
        print >> out, obs, p, maxpost
    o_stream.close()
    return out.getvalue()


def main(args):
    filenames = (args.out_forward, args.out_scaling, args.out_backward)
    # aggregate and validate the model parameters
    model = DGRP.Model()
    model.from_fieldstorage(args)
    # see how the states interact with the observations
    states = (
            model.get_recent_state(),
            model.get_ancient_state(),
            model.get_misaligned_state(args.misalignment_effect),
            model.get_garbage_state())
    # define the transition object
    nstates = len(states)
    prandom = min(1.0, (nstates / (nstates - 1.0)) / args.region_size)
    T = TransitionMatrix.UniformTransitionObject(prandom, nstates)
    # make the hmm
    hmm = ExternalHMM.ExternalModel(T, states, filenames)
    converter = lineario.IntTupleConverter()
    o_stream = lineario.SequentialDiskIO(converter, args.obsfile)
    hmm.init_dp(o_stream)
    o_stream.open_read()
    for p, obs in itertools.izip(hmm.posterior(), o_stream.read_forward()):
        p_recent, p_ancient, p_misaligned, p_garbage = p
        # get the prior probability of polymorphism conditional on state
        p_recent_AA = states[0].get_posterior_distribution(obs)[2]
        p_ancient_AA = states[1].get_posterior_distribution(obs)[2]
        # compute the posterior probability of a polymorphism
        posterior_polymorphism = 0
        posterior_polymorphism += p_recent * p_recent_AA
        posterior_polymorphism += p_ancient * p_ancient_AA
        # Given that a polymorphism occurred,
        # get the probability distribution over the
        # three non-reference nucleotides.
        r = model.seqerr
        log_Pr = math.log(r/4.0)
        log_PA = math.log(1 - 3*r/4.0)
        logs = [
                obs[1]*log_PA + obs[2]*log_Pr + obs[3]*log_Pr,
                obs[1]*log_Pr + obs[2]*log_PA + obs[3]*log_Pr,
                obs[1]*log_Pr + obs[2]*log_Pr + obs[3]*log_PA]
        condmaxpost = math.exp(max(logs) - scipy.misc.logsumexp(logs))
        # get the posterior probability distribution
        maxpost = posterior_polymorphism * condmaxpost
        # show the annotation for this position
        annotation = list(obs) + list(p) + [maxpost]
        print '\t'.join(str(x) for x in annotation)
    o_stream.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--out_forward',
            help='forward vectors are written to this file')
    parser.add_argument('--out_backward',
            help='backward vectors are written to this file')
    parser.add_argument('--out_scaling',
            help='scaling factors are written to this file')
    parser.add_argument('--misalignment_effect', type=float, default=0.5,
            help='misalignment branch length')
    parser.add_argument('--x', type=float, default=0.1,
            help='reference branch length')
    parser.add_argument('--y', type=float, default=0.01,
            help='line branch length')
    parser.add_argument('--z', type=float, default=0.0001,
            help='mutant branch length')
    parser.add_argument('--low', type=int, default=2,
            help='low random coverage per base')
    parser.add_argument('--med', type=int, default=20,
            help='medium random coverage per base')
    parser.add_argument('--high', type=int, default=1000,
            help='high random coverage per base')
    parser.add_argument('--seqerr', type=float, default=0.1,
            help='sequencing error')
    parser.add_argument('--nomcoverage', type=int, default=20,
            help='nominal total coverage per position')
    parser.add_argument('--kmulticoverages', type=int, default=4,
            help='allowed multiples of nominal coverage')
    parser.add_argument('--region_size', type=int, default=1000,
            help='expected contiguous region lengths')
    parser.add_argument('obsfile')
    args = parser.parse_args()
    main(args)
