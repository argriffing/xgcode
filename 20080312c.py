"""Sample a codon alignment given a tree and a Direct Protein mixture model.

The branch length is the expected number of codon substitutions on a branch.
See "Population Genetics Without Intraspecific Data" by Thorne et al.
for more information about the Direct Protein model.
"""

from SnippetUtil import HandlingError
import Newick
import PhyLikelihood
import DirectProtein
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree',
                formatted_tree_string),
            Form.MultiLine('model', 'Direct Protein mixture model',
                DirectProtein.get_sample_xml_string().strip()),
            Form.Integer('ncols', 'sample this many codon columns',
                100, low=1, high=1000)]
    return form_objects

def get_form_out():
    return FormOut.CodonFasta()

def get_response_content(fs):
    # get the newick string
    try:
        tree = Newick.parse(fs.tree, Newick.NewickTree)
        tree.assert_valid()
    except Newick.NewickSyntaxError as e:
        raise HandlingError(e)
    # get the normalized Direct RNA mixture model
    mixture_model = DirectProtein.deserialize_mixture_model(fs.model)
    mixture_model.normalize()
    # simulate the alignment
    try:
        alignment = PhyLikelihood.simulate_alignment(tree,
                mixture_model, fs.ncols)
    except PhyLikelihood.SimulationError as e:
        raise HandlingError(e)
    # get the alignment
    arr = []
    for node in tree.gen_tips():
        arr.append(alignment.get_fasta_sequence(node.name))
    # return the alignment string
    return '\n'.join(arr) + '\n'
