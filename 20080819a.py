"""Sample an aligned pair of nucleotide sequences using the F84 model.

The evolutionary distance is the expected number of substitutions per site.
The F84 evolutionary model is defined in the paper
"Maximum Likelihood Phylogenetic Estimation from DNA Sequences
with Variable Rates over Sites: Approximate Methods"
by Ziheng Yang in J Mol Evol 1994.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import F84
import PairLikelihood
import Fasta
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('length', 'sequence length',
                60, low=2, high=1000),
            Form.Float('distance', 'evolutionary distance',
                1.0, low_exclusive=0),
            Form.Float('kappa', 'kappa', 2.0, low_inclusive=0),
            Form.Float('A', 'weight of A', 1.0, low_inclusive=0),
            Form.Float('C', 'weight of C', 1.0, low_inclusive=0),
            Form.Float('G', 'weight of G', 1.0, low_inclusive=0),
            Form.Float('T', 'weight of T', 1.0, low_inclusive=0)]
    return form_objects

def get_form_out():
    return FormOut.NucleotideFasta()

def get_response_content(fs):
    # read the nucleotide weights
    nt_weights = [fs.A, fs.C, fs.G, fs.T]
    # convert the nucleotide weights to probabilities
    nt_probs = [x / float(sum(nt_weights)) for x in nt_weights]
    # Assert that the kappa value and the nucleotide
    # probabilities are compatible.
    A, C, G, T = nt_probs
    R = float(A + G)
    Y = float(C + T)
    if R <= 0:
        raise HandlingError('the frequency of a purine must be positive')
    if Y <= 0:
        raise HandlingError('the frequency of a pyrimidine must be positive')
    if fs.kappa <= max(-Y, -R):
        msg_a = 'kappa must be greater than max(-R, -Y) '
        msg_b = 'where R and Y are the purine and pyrimidine frequencies'
        raise HandlingError(msg_a + msg_b)
    # Create the rate matrix object
    # which is automatically scaled to a rate of 1.0.
    model = F84.create_rate_matrix(fs.kappa, nt_probs)
    # simulate a pair of sequences
    sequence_pair = PairLikelihood.simulate_sequence_pair(
            fs.distance, model, fs.length)
    # convert the pair of sequences to an alignment object
    aln = StringIO()
    print >> aln, '>first'
    print >> aln, ''.join(sequence_pair[0])
    print >> aln, '>second'
    print >> aln, ''.join(sequence_pair[1])
    return Fasta.Alignment(StringIO(aln.getvalue())).to_fasta_string() + '\n'
