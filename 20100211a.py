"""Analyze a single nucleotide site of a pileup.

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

from SnippetUtil import HandlingError
import Form
import DGRP
import ReadCoverageRef
import Util

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

def get_form():
    """
    @return: the body of a form
    """
    default_counts = [
            ('A', '20'),
            ('C', '20'),
            ('G', '3'),
            ('T', '1')]
    form_objects = [
            Form.MultiLine('param_field', 'parameters',
                '\n'.join('\t'.join(p) for p in g_default_params)),
            Form.MultiLine('count_field', 'nucleotide read counts',
                '\n'.join('\t'.join(p) for p in default_counts)),
            Form.SingleLine('ref', 'reference nucleotide', 'A')]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    # define the observation
    lines = Util.get_stripped_lines(StringIO(fs.count_field))
    if len(lines) != 4:
        raise Exception('expected 4 nucleotide lines')
    rows = [line.split() for line in lines]
    if any(len(r)!=2 for r in rows):
        raise Exception('nucleotide read count syntax error')
    typed_rows = [(nt, int(count)) for nt, count in rows]
    nt_to_count = dict(typed_rows)
    if set(nt_to_count.keys()) != set('ACGT'):
        raise Exception('invalid nucleotide in the count field')
    ref = fs.ref.strip()
    if ref not in list('ACGT'):
        raise Exception('invalid reference nucleotide')
    obs = [nt_to_count[ref]]
    obs += [nt_to_count[nt] for nt in 'ACGT' if nt != ref]
    # read the parameters
    lines = Util.get_stripped_lines(StringIO(fs.param_field))
    rows = [line.split() for line in lines]
    if any(len(r)!=2 for r in rows):
        raise Exception('parameter syntax error')
    default_params, default_values = zip(*g_default_params)
    param_to_value = dict(rows)
    if set(default_params) != set(param_to_value.keys()):
        raise Exception('invalid or missing parameter')
    x = float(param_to_value['x'])
    y = float(param_to_value['y'])
    z = float(param_to_value['z'])
    low = int(param_to_value['low'])
    med = int(param_to_value['med'])
    high = int(param_to_value['high'])
    seqerr = float(param_to_value['seqerr'])
    nom = int(param_to_value['nomcoverage'])
    kmulti = int(param_to_value['kmulticoverages'])
    # create the three states
    recent = ReadCoverageRef.HMMRecent(x, y, z, seqerr, nom, kmulti)
    ancient = ReadCoverageRef.HMMAncient(x, y, z, seqerr, nom, kmulti)
    garbage = ReadCoverageRef.HMMGarbage(low, med, high)
    # see how the three states interact with the observation
    states = (recent, ancient, garbage)
    names = ('recent', 'ancient', 'garbage')
    likelihoods = [s.get_likelihood(obs) for s in states]
    log_likelihoods = [s.get_log_likelihood(obs) for s in states]
    # report the log likelihoods
    print >> out, 'log likelihoods:'
    for name, loglik in zip(names, log_likelihoods):
        print >> out, '%s:\t%s' % (name, loglik)
    print >> out
    # report the posterior distribution
    posterior = [v/sum(likelihoods) for v in likelihoods]
    print >> out, 'posterior distribution:'
    for name, p in zip(names, posterior):
        print >> out, '%s:\t%s' % (name, p)
    print >> out
    return [('Content-Type', 'text/plain')], out.getvalue().strip()
