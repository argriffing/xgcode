"""Given a newick tree, use JC69 to sample aligned nt sequences at the leaves.

Given a newick tree, use JC69 to sample aligned
nucleotide sequences at the leaves.
JC69 is the simple continuous time Markov model
proposed by Jukes and Cantor in 1969.
The sequence order field may be left empty
if the order of the sequences in the FASTA output is unimportant.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Newick
import Util
import Fasta
import JC69
import Form
import FormOut

#FIXME use const data

# This tree is in an arXiv paper by Lior Pachter.
g_newick = '(((a:0.05, b:0.05):0.15, c:0.2):0.8, x:1.0, (((m:0.05, n:0.05):0.15, p:0.2):0.8, y:1.0):1.0);'

def get_form():
    """
    @return: the body of a form
    """
    # define the newick string and the ordered labels
    tree_string = g_newick
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    labels = list('xyabcmnp')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('order', 'sequence order', '\n'.join(labels)),
            Form.Integer('length', 'sequence length', 60, low=2)]
    return form_objects

def get_form_out():
    return FormOut.NucleotideFasta()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    # get the sequence order if it exists
    ordered_names = Util.get_stripped_lines(StringIO(fs.order))
    if ordered_names:
        observed_name_set = set(ordered_names)
        expected_name_set = set(node.get_name() for node in tree.gen_tips())
        extra_names = observed_name_set - expected_name_set
        missing_names = expected_name_set - observed_name_set
        if extra_names:
            msg_a = 'the list of ordered names includes these names '
            msg_b = 'not found in the tree: %s' % str(tuple(extra_names))
            raise HandlingError(msg_a + msg_b)
        if missing_names:
            msg_a = 'the tree includes these names not found in the list '
            msg_b = 'of ordered names: %s' % str(tuple(missing_names))
            raise HandlingError(msg_a + msg_b)
    else:
        ordered_names = list(tip.get_name() for name in tree.gen_tips())
    # do the sampling
    sampled_sequences = JC69.sample_sequences(tree, ordered_names, fs.length)
    alignment = Fasta.create_alignment(ordered_names, sampled_sequences)
    # begin the response
    out = StringIO()
    print >> out, alignment.to_fasta_string()
    # return the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
