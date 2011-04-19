"""
Read Newick data.
This module reads and writes Newick data
in a way that does not depend much on the underlying tree representation.
It is closely related to the Newick module
which imports many of these functions.
"""

import unittest

# This is a rooted version of the tree on page 573 of Inferring Phylogenies by Joseph Felsenstein.
rooted_example_tree = "(((((((A:4, B:4):6.125, C:5.1):8, D:6):3, E:21):10, ((F:4, G:12):14, H:8):13):13, ((I:5, J:2):30, (K:11, L:11):2):17):4, M:56);"

# This tree is on page 573 of Inferring Phylogenies by Joseph Felsenstein.
daylight_example_tree = "((((((A:4, B:4):6, C:5):8, D:6):3, E:21):10, ((F:4, G:12):14, H:8):13):13, ((I:5, J:2):30, (K:11, L:11):2):17, M:56);"

# This tree is provided by paml.
brown_example_tree = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'

# This is from Ben Redelings, and probably indirectly from somewhere else.
ben_example_tree = """
(Dif_Fly_[d.melanogaster]:0.824935,(dorsal-mosquito.aedes1:0.340352,horseshoe-crab:0.347596,(Dorsal_Fly_[d.melanogaster]:0.003629,dorsal-d.sim:0.006374):0.207850,(RelA-beetle:0.342349,dorsal-honeybee2:0.149771):0.088256):0.190091,(squid:0.497853,(Urchin:0.989007,(((Rel_Sea_Peach_[h.roretzi]:0.283520,(Rel_Sea_Squirt_[c.intestinalis]:0.069185,Rel_Sea_Squirt_[c.savignyi]:0.085474):0.339533):0.261828,(Rel_Lancelet_[b.florida]:0.483145,(RelB_Zebrafish:0.406528,(RelB_Frog:0.345967,(RelB_Mouse:0.130419,RelB_chicken:0.194183):0.142151):0.225204):0.263858):0.090565,(Rel_Sea_Lamprey_[p.marinus]:0.350536,((RelA_Zebrafish:0.229731,(RelA_Frog:0.118917,(RelA_Mouse:0.125854,RelA_chicken:0.312344):0.073130):0.089640):0.171363,(cRel_Zebrafish:0.181109,(cRel_Mouse:0.148298,(cRel_Frog:0.271300,cRel_chicken:0.108538):0.023230):0.066827):0.104448):0.177601):0.181920):0.111531,(Relish_Fly_[d.melanogaster]:1.211823,(Urchin2:0.359060,((NFkB_Sea_Squirt_[c.intestinalis]:0.026731,NFkB_Sea_Squirt_[c.savignyi]:0.028602):0.405836,(p50_Frog:0.120142,(p50_Mouse:0.047345,p50_chicken:0.033400):0.079139):0.228715,(p52_Zebrafish:0.231553,(p52_Frog:0.142828,(p52_Mouse:0.116143,p52_chicken:0.067973):0.059010):0.128246):0.208689):0.102472):0.153007):0.196171):0.095239):0.099623):0.183129);
"""

# This is from the interactive tree of life.
itol_archaea_tree = """
(((Pyrobaculum_aerophilum:0.55514,(Aeropyrum_pernix:0.43380,(Sulfolobus_tokodaii:0.17389,Sulfolobus_solfataricus:0.18962)100:0.33720)94:0.09462)100:0.12018,Nanoarchaeum_equitans:0.81078)100:0.15444,((Thermoplasma_acidophilum:0.09785,Thermoplasma_volcanium:0.10412)100:0.66151,(((Pyrococcus_furiosus:0.02366,(Pyrococcus_abyssi:0.02179,Pyrococcus_horikoshii:0.02786)100:0.02239)100:0.36220,((Methanopyrus_kandleri:0.35331,Methanobacterium_thermautotrophicum:0.36583)99:0.07446,(Methanococcus_jannaschii:0.13226,Methanococcus_maripaludis:0.28592)100:0.23828)100:0.06284)51:0.04469,(Archaeoglobus_fulgidus:0.34660,((Methanosarcina_acetivorans:0.02602,Methanosarcina_mazei:0.03087)100:0.30588,Halobacterium_sp._NRC-1:0.61597)100:0.12801)100:0.10395)62:0.06815)99:0.11833)Archaea;
"""

class NewickSyntaxError(Exception):
    pass

def _lex_newick_raw(s):
    """
    Split a newick string into punctuation and non-punctuation.
    """
    punctuation = tuple('():,;[]') 
    symbol = ''
    for c in s:
        if c in punctuation:
            yield symbol
            yield c
            symbol = ''
        else:
            symbol += c
    yield symbol

def _lex_newick(s):
    """
    Split a newick string.
    Split a newick string into meaningful symbols
    separated by punctuation / non-punctuation boundaries.
    Note that this method allows whitespace inside node names
    but not leading or terminal whitespace.
    Also this function removes bracketed symbols.
    """
    bracket_depth = 0
    for symbol in _lex_newick_raw(s):
        symbol = symbol.strip()
        if not symbol:
            continue
        if symbol == '[':
            bracket_depth += 1
        elif symbol == ']':
            bracket_depth -= 1
        elif bracket_depth == 0:
            yield(symbol)
        if bracket_depth < 0:
            raise NewickSyntaxError('unmatched bracket')
    if bracket_depth != 0:
        raise NewickSyntaxError('unbalanced brackets')


# In this experimental section we use a new API.
# The new API is flatter and expects a tree object (not a factory)
# from the caller, where the tree object has
# four member functions:  v=create_root(), v=create_child(parent),
# set_branch_length(v,length), and set_root(v) .
# A compatibility level will be added.

def _pnh_blen(tree, symbols, index, node):
    """
    This uses the simplified API.
    Parse Newick Helper.
    This is called when we are at a ":" symbol.
    @param tree: the object provided by the user
    @param node: a node in the user tree
    @return: next index
    """
    blen_str = symbols[index+1]
    try:
        branch_length = float(blen_str)
    except ValueError:
        msg = 'found "%s" instead of a branch length' % blen_str
        raise NewickSyntaxError(msg)
    tree.set_branch_length(node, branch_length)
    return index+2

def _pnh(tree, symbols, index):
    """
    This uses the simplified API.
    Parse Newick Helper.
    @param tree: the object provided by the user
    @param symbols: the output of the newick lexer
    @param index: the current position in the symbol list
    @return: (root node, next index)
    """
    if symbols[-1] != ';':
        msg = 'the newick symbol list should end with a semicolon'
        raise NewickSyntaxError(msg)
    if index >= len(symbols):
        raise NewickSyntaxError('premature string termination')
    root = tree.create_root()
    if symbols[index] == ',':
        # Hitting this symbol when we expect a new node means that
        # an unnamed node with no branch length was found.
        return (root, index)
    elif symbols[index] == ')':
        # Hitting this symbol when we expect a new node means that
        # an unnamed node with no branch length was found
        # and that the list terminates.
        return (root, index)
    elif symbols[index] == ':':
        # Hitting this symbol when we expect a new node means that
        # an unnamed node with a branch length was found.
        next_index = _pnh_blen(tree, symbols, index, root)
        return (root, next_index)
    elif symbols[index] == '(':
        # Hitting this symbol when we expect a new node means that
        # one or more child nodes must be created
        # in addition to the root node.
        index += 1
        while True:
            child, index = _pnh(tree, symbols, index)
            tree.add_child(root, child)
            if symbols[index] == ',':
                index += 1
            elif symbols[index] == ')':
                if symbols[index+1] not in list(':;(),'):
                    root.add_name(symbols[index+1])
                    index += 1
                if symbols[index+1] == ':':
                    next_index = _pnh_blen(tree, symbols, index+1, root)
                else:
                    next_index = index+1
                return (root, next_index)
            else:
                msg_a = 'found "%s" instead of ' % symbols[index]
                msg_b = 'a comma or a closing parenthesis'
                raise NewickSyntaxError(msg_a + msg_b)
    elif symbols[index] == ';':
        raise NewickSyntaxError('found the ";" terminator prematurely')
    else:
        # Hitting an unrecognized symbol means that
        # we are starting a named node.
        tree.set_name(root, symbols[index])
        if symbols[index+1] == ':':
            next_index = _pnh_blen(tree, symbols, index+1, root)
        else:
            next_index = index+1
        return (root, next_index)

def parse_simplified(s, tree):
    """
    This is a new parse function.
    Note that this takes a tree instead of a tree factory.
    @param s: the newick tree string to be parsed
    @param tree: a newly created tree
    @return: the tree
    """
    if not s:
        raise NewickSyntaxError('empty tree string')
    symbols = list(_lex_newick(s))
    if symbols.count('(') != symbols.count(')'):
        raise NewickSyntaxError('parenthesis mismatch')
    if not symbols[-1] == ';':
        msg_a = 'the newick symbol list should end with a semicolon: '
        msg_b = str(symbols)
        raise NewickSyntaxError(msg_a + msg_b)
    root, index = _pnh(tree, symbols, 0)
    if index >= len(symbols):
        msg = 'the parser tried to use too much of the newick string'
        raise NewickSyntaxError(msg)
    if index < len(symbols) - 1:
        msg = 'the parser did not use the whole newick string'
        raise NewickSyntaxError(msg)
    tree.set_root(root)
    return tree


# In this section we add API wrappers for the old interface.

class _T_wrapper:
    def __init__(self, tree):
        self.tree = tree
    def create_root(self):
        return self.tree.NodeFactory()
    def add_child(self, parent, child):
        parent.add_child(child)
        child.set_parent(parent)
    def set_branch_length(self, node, blen):
        node.set_branch_length(blen)
    def set_name(self, node, name):
        node.add_name(name)
    def set_root(self, node):
        self.tree.set_root(node)

def parse(s, tree_factory):
    """
    This should act like the old parse function.
    @param s: the newick tree string to be parsed
    @param tree_factory: a callable that returns a tree given a root
    @return: the tree
    """
    wrapped_tree = parse_simplified(s, _T_wrapper(tree_factory()))
    return wrapped_tree.tree


# Remaining functions.

def _get_name_string(node):
    """
    Get a string representation of the node.
    This is the name of the node if one exists,
    but otherwise it a string representation of the node id.
    @return: the string representation of the node
    """
    name = node.get_name()
    if name is None:
        return '<node %d>' % id(node)
    else:
        return name

def _get_blen_string(node):
    """
    To save space trim each trailing zero.
    """
    branch_length = node.get_branch_length()
    if branch_length is None:
        return None
    if int(branch_length) == branch_length:
        return str(int(branch_length))
    return str(branch_length).rstrip('0')


def _get_multiline_newick_string_helper(node, nlevels, depth):
    """
    Yield lines.
    @param node: a newick node
    @param nlevels: the number of levels to expand
    """
    subtree_lists = []
    for child in node.gen_children():
        lines = list(_get_multiline_newick_string_helper(
            child, nlevels, depth+1))
        subtree_lists.append(lines)
    if depth < nlevels:
        # yield multiple lines with indentation
        if subtree_lists:
            yield '('
            for i, subtree_list in enumerate(subtree_lists):
                for j, line in enumerate(subtree_list):
                    # indent every line
                    outline = '  ' + line
                    # put a comma at the end of some of the lines
                    jfinal = (j == len(subtree_list) - 1)
                    ifinal = (i == len(subtree_lists) - 1)
                    if jfinal and not ifinal:
                        outline += ','
                    yield outline
            internal_name_string = ''
            if node.name is not None:
                internal_name_string = node.name
            base_string = ')' + internal_name_string
        else:
            base_string = _get_name_string(node)
        if node.get_branch_length() is not None:
            blen_string = _get_blen_string(node)
            yield '%s:%s' % (base_string, blen_string)
        else:
            yield base_string
    else:
        # yield a single line with no indentation
        if subtree_lists:
            subtree_strings = [lines[0] for lines in subtree_lists]
            internal_name_string = ''
            if node.name is not None:
                internal_name_string = node.name
            base_string = ''.join([
                '(', ', '.join(subtree_strings), ')', internal_name_string])
        else:
            base_string = _get_name_string(node)
        blen_string = _get_blen_string(node)
        if node.get_branch_length() is not None:
            yield '%s:%s' % (base_string, blen_string)
        else:
            yield base_string

def get_newick_string(tree):
    """
    @param tree: a newick tree
    @return: a single line newick string
    """
    return get_multiline_newick_string(tree, 0)

def get_multiline_newick_string(tree, nlevels):
    """
    @param tree: a newick tree
    @param nlevels: the number of levels to expand
    @return: a multi-line newick string
    """
    lines = list(_get_multiline_newick_string_helper(tree.root, nlevels, 0))
    return '\n'.join(lines) + ';'

def get_narrow_newick_string(tree, maxwidth):
    """
    @param tree: a newick tree
    @param maxwidth: the maximum allowed line length including indentation
    @return: a multi-line newick string
    """
    for nlevels in range(tree.get_height() + 1):
        lines = list(_get_multiline_newick_string_helper(tree.root, nlevels, 0))
        if max(len(line) for line in lines) <= maxwidth:
            break
    return '\n'.join(lines) + ';'





class TestNewick(unittest.TestCase):

    # override TreeFactory in derived classes
    TreeFactory = None

    def assert_correct_syntax(self, s):
        if not self.TreeFactory:
            return
        try:
            parse(s, self.TreeFactory)
        except NewickSyntaxError, e:
            self.fail(e)

    # The following four tests assert the correctness of strings provided by Felsenstein.
    # http://evolution.genetics.washington.edu/phylip/newicktree.html

    def test_multifurcation(self):
        self.assert_correct_syntax('(B,(A,C,E),D);')

    def test_degenerate_nodes(self):
        self.assert_correct_syntax('(,(,,),);')

    def test_named_ancestor(self):
        self.assert_correct_syntax('(B:6.0,(A:5.0,C:3.0,E:4.0)Ancestor1:5.0,D:11.0);')

    def test_itol_archaea(self):
        self.assert_correct_syntax(itol_archaea_tree)

    def test_felsenstein_strings(self):
        valid_newick_strings = (
                '((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154);',
                '(Bovine:0.69395,(Gibbon:0.36079,(Orang:0.33636,(Gorilla:0.17147,(Chimp:0.19268, Human:0.11927):0.08386):0.06124):0.15057):0.54939,Mouse:1.21460):0.10;',
                '(Bovine:0.69395,(Hylobates:0.36079,(Pongo:0.33636,(G._Gorilla:0.17147, (P._paniscus:0.19268,H._sapiens:0.11927):0.08386):0.06124):0.15057):0.54939, Rodent:1.21460);',
                'A;',
                '((A,B),(C,D));',
                '(Alpha,Beta,Gamma,Delta,,Epsilon,,,);'
                )
        for s in valid_newick_strings:
            self.assert_correct_syntax(s)

    def test_get_multiline_newick_string_shallow(self):
        if not self.TreeFactory:
            return
        tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
        tree = parse(tree_string, self.TreeFactory)
        observed = get_multiline_newick_string(tree, 1)
        arr = [
                '(',
                '  ((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7,',
                '  Orangutan:0.4,',
                '  Gibbon:0.5',
                ');']
        expected = '\n'.join(arr)
        self.assertEquals(observed, expected)

    def test_get_multiline_newick_string_deep(self):
        if not self.TreeFactory:
            return
        tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
        tree = parse(tree_string, self.TreeFactory)
        observed = get_multiline_newick_string(tree, 2)
        arr = [
                '(',
                '  (',
                '    (Human:0.1, Chimpanzee:0.2):0.8,',
                '    Gorilla:0.3',
                '  ):0.7,',
                '  Orangutan:0.4,',
                '  Gibbon:0.5',
                ');']
        expected = '\n'.join(arr)
        self.assertEquals(observed, expected)

    def test_get_narrow_newick_string(self):
        if not self.TreeFactory:
            return
        tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
        tree = parse(tree_string, self.TreeFactory)
        observed = get_narrow_newick_string(tree, 60)
        arr = [
                '(',
                '  ((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7,',
                '  Orangutan:0.4,',
                '  Gibbon:0.5',
                ');']
        expected = '\n'.join(arr)
        self.assertEquals(observed, expected)

    def test_get_narrow_newick_string_long(self):
        if not self.TreeFactory:
            return
        tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
        tree = parse(tree_string, self.TreeFactory)
        observed = get_narrow_newick_string(tree, 1000)
        self.assertEquals(tree_string, observed)

    def test_get_narrow_newick_string_short(self):
        if not self.TreeFactory:
            return
        tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
        tree = parse(tree_string, self.TreeFactory)
        observed_a = get_narrow_newick_string(tree, 2)
        observed_b = get_narrow_newick_string(tree, 5)
        self.assertEquals(observed_a, observed_b)

    def test_printing_without_branch_lengths(self):
        if not self.TreeFactory:
            return
        tree_string = '(((a, b), c), x, (((m, n), p), y));'
        tree = parse(tree_string, self.TreeFactory)
        observed = get_newick_string(tree)
        expected = tree_string
        self.assertEquals(observed, expected)


if __name__ == '__main__':
    """
    This module has only generic functions.
    These tests should be run from modules with concrete implementations.
    """
    pass
