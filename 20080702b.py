"""Evaluate a bipartition function given a set of trees.

The exact bipartition criterion is a matrix function by Eric Stone.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Util
import FelTree
import NewickIO
import Clustering
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define some default tree lines
    default_tree_lines = [
            '(a:1.062, c:0.190, (d:1.080, b:2.30):2.112);',
            '((d:0.083, b:0.614):0.150, e:0.581, (c:1.290, a:0.359):1.070);',
            '((b:1.749, d:0.523):0.107, e:1.703, (a:0.746, c:0.070):4.025);']
    # define the form objects
    form_objects = [
            Form.MultiLine('trees', 'newick trees (one tree per line)', '\n'.join(default_tree_lines)),
            Form.RadioGroup('criterion', 'bipartition function', [
                Form.RadioItem('exact', 'exact criterion', True),
                Form.RadioItem('sign', 'spectral sign approximation'),
                Form.RadioItem('threshold', 'spectral threshold approximation'),
                Form.RadioItem('nj', 'neighbor joining criterion'),
                Form.RadioItem('random', 'random bipartition')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # get the newick trees.
    trees = []
    for tree_string in Util.stripped_lines(StringIO(fs.trees)):
        # parse each tree and make sure that it conforms to various requirements
        tree = NewickIO.parse(tree_string, FelTree.NewickTree)
        tip_names = [tip.get_name() for tip in tree.gen_tips()]
        if len(tip_names) < 4:
            raise HandlingError('expected at least four tips but found ' + str(len(tip_names)))
        if any(name is None for name in tip_names):
            raise HandlingError('each terminal node must be labeled')
        if len(set(tip_names)) != len(tip_names):
            raise HandlingError('each terminal node label must be unique')
        trees.append(tree)
    # read the criterion string, creating the splitter object
    if fs.exact:
        splitter = Clustering.StoneExactDMS()
    elif fs.sign:
        splitter = Clustering.StoneSpectralSignDMS()
    elif fs.threshold:
        splitter = Clustering.StoneSpectralThresholdDMS()
    elif fs.nj:
        splitter = Clustering.NeighborJoiningDMS()
    elif fs.random:
        splitter = Clustering.RandomDMS()
    # assert that the computation is fast
    complexity = 0
    for tree in trees:
        n = len(list(tree.gen_tips()))
        complexity += n * splitter.get_complexity(n)
    if complexity > 1000000:
        raise HandlingError('this computation would take too long')
    # evaluate the bipartition of each tree based on its distance matrix
    informative_split_count = 0
    degenerate_split_count = 0
    invalid_split_count = 0
    for tree in trees:
        tips = list(tree.gen_tips())
        n = len(tips)
        distance_matrix = tree.get_distance_matrix()
        index_selection = splitter.get_selection(distance_matrix)
        tip_selection = [tips[i] for i in index_selection]
        n_selection = len(tip_selection)
        n_complement = n - n_selection
        if min(n_selection, n_complement) < 2:
            degenerate_split_count += 1
        else:
            if tree.get_split_branch(tip_selection):
                informative_split_count += 1
            else:
                invalid_split_count += 1
    # define the response
    out = StringIO()
    print >> out, informative_split_count, 'informative splits'
    print >> out, degenerate_split_count, 'degenerate splits'
    print >> out, invalid_split_count, 'invalid splits'
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
