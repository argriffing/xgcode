"""Analyze multiple sequential nucleotide sites of a pileup.

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

from SnippetUtil import HandlingError
import Form
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
            Form.Float('misalignment_effect', 'misalignment_effect',
                '0.5', low_inclusive=0),
            Form.MultiLine('param_field', 'parameters',
                '\n'.join('\t'.join(p) for p in g_default_params)),
            Form.MultiLine('data_field', 'data',
                '\n'.join(data_lines))]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # define the three HMM states
    lines = Util.get_stripped_lines(StringIO(fs.param_field))
    model = DGRP.Model()
    model.from_lines(lines)
    # see how the three states interact with the observations
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
    obs_a = lineario.SequentialStringIO(converter, fs.data_field)
    obs_b = lineario.SequentialStringIO(converter, fs.data_field)
    obs_b.open_read()
    hmm.init_dp(obs_a)
    state_header = ('recent', 'ancient', 'misaligned', 'garbage')
    obs_header = ('ref', 'nonref-x', 'nonref-y', 'nonref-z')
    print >> out, obs_header, state_header
    for p, obs in itertools.izip(hmm.posterior(), obs_b.read_forward()):
        print >> out, obs, p
    obs_b.close()
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
