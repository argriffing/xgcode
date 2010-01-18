"""Given a newick tree, calculate weights for a subset of the tips.
"""

import StringIO

from SnippetUtil import HandlingError
import Newick
import LeafWeights
import Util
import Form

# This tree is from UCSC.
g_ucsc_tree_string = """
(((((((((
(human_hg18:0.00669,chimp_panTro1:0.00757):0.0243,
macaque_rheMac2:0.0592):0.0240,
((rat_rn4:0.0817,mouse_mm8:0.0770):0.229,
rabbit_oryCun1:0.207):0.107):0.0230,
(cow_bosTau2:0.159,dog_canFam2:0.148):0.0394):0.0285,
armadillo_dasNov1:0.150):0.0160,
(elephant_loxAfr1:0.105,tenrec_echTel1:0.260):0.0404):0.218,
monodelphis_monDom4:0.371):0.189,
chicken_galGal2:0.455):0.123,
xenopus_xenTro1:0.782):0.156,
((tetraodon_tetNig1:0.199,fugu_fr1:0.239):0.493,
zebrafish_danRer3:0.783):0.156);
"""

def get_form():
    """
    @return: the body of a form
    """
    # define the default tree string
    tree_string = g_ucsc_tree_string
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # default name selection
    default_name_selection = ('human_hg18', 'chimp_panTro1', 'rabbit_oryCun1')
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('selection', 'selected tip names', '\n'.join(default_name_selection)),
            Form.RadioGroup('method', 'weighting method', [
                Form.RadioItem('stone', 'use the method of Stone and Sidow', True),
                Form.RadioItem('thompson', 'use the method of Thompson et al.')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the tree
    tree = Newick.parse(fs.tree, Newick.NewickTree)
    tree.assert_valid()
    tree.add_branch_lengths()
    if tree.has_negative_branch_lengths():
        raise HandlingError('calculating weights for a tree with negative branch lengths is not implemented')
    # get the selected names
    selection = list(Util.stripped_lines(StringIO.StringIO(fs.selection)))
    selected_name_set = set(selection)
    possible_name_set = set(node.get_name() for node in tree.gen_tips())
    extra_names = selected_name_set - possible_name_set
    if extra_names:
        raise HandlingError('the following selected names are not valid tips: %s' % str(tuple(extra_names)))
    # prune the tree 
    for name in set(node.name for node in tree.gen_tips()) - set(selection): 
        try: 
            node = tree.get_unique_node(name) 
        except NewickSearchError, e: 
            raise HandlingError(e) 
        tree.prune(node)
    # get the weights
    if fs.stone:
        name_weight_pairs = LeafWeights.get_stone_weights(tree)
    elif fs.thompson:
        name_weight_pairs = LeafWeights.get_thompson_weights(tree)
    # report the weights
    out = StringIO.StringIO()
    for name, weight in name_weight_pairs:
        print >> out, '%s: %f' % (name, weight)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
