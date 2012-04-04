"""
Sample some data so that BEAST can later guess the rate autocorrelation.[UNFINISHED]

BEAST defines this as the covariance between parent lineage rates
and child lineage rates.
In its covariance calculation, BEAST uses Bessel's correction
referred to in its comments as a correction so that the estimate is ML,
but this is a mistake because the correction is actually so that
the estimate is unbiased.
"""

from StringIO import StringIO
import string
import math
import random
from collections import defaultdict

import Form
import FormOut
import kingman
import FtreeIO
import Ftree

#TODO use lxml

g_tree_model_defn = """
<treeModel id="treeModel1">
    <newick idref="startingTree"/>
    <rootHeight>
        <parameter id="treeModel1.rootHeight"/>
    </rootHeight>
    <nodeHeights internalNodes="true">
        <parameter id="treeModel1.internalNodeHeights"/>
    </nodeHeights>
    <nodeHeights internalNodes="true" rootNode="true">
        <parameter id="treeModel1.allInternalNodeHeights"/>
    </nodeHeights>
</treeModel>
""".strip()

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.Integer('nleaves', 'sample this many leaf taxa',
                10, low=3, high=26),
            Form.Integer('nsites', 'sample this many alignment columns',
                200, low=1, high=2000),
            Form.Float('epsrate', 'larger values mean more rate variation',
                0.02, low_exclusive=0, high_exclusive=1),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.BeastXml('simulation.xml')

def sample_jc_column(R, B):
    """
    Sample a column of a Jukes-Cantor alignment.
    @param R: Ftree directed topology
    @param B: branch lengths in expected number of substitutions
    @return: a map from vertex to nucleotide
    """
    acgt = 'ACGT'
    v_to_nt = {}
    v_to_source = R_to_v_to_source(R)
    for v in R_to_preorder(R):
        p = v_to_source.get(v, None)
        if p is None:
            v_to_nt[v] = random.choice(acgt)
        else:
            d = B[frozenset([v, p])]
            p_randomize = 1.0 - math.exp(-(4.0 / 3.0) * d)
            if random.random() < p_randomize:
                v_to_nt[v] = random.choice(acgt)
            else:
                v_to_nt[v] = v_to_nt[p]
    return v_to_nt

def R_to_v_to_leaves(R):
    """
    @return: map from v to its set of subtree leaves
    """
    v_to_leaves = defaultdict(set)
    for path in kingman.get_paths_to_root(R):
        for v in path:
            v_to_leaves[v].add(path[0])
    return v_to_leaves

def get_mrca_subset_defn(N, v, leaves):
    """
    Get the first chunk of http://beast.bio.ed.ac.uk/Tutorial_3.1
    @param N: map from vertex to name
    @param v: the internal vertex
    @param leaves: the leaves below v
    """
    out = StringIO()
    print >> out, '<taxa id="taxa_subset_%s">' % N[v]
    for leaf in sorted(leaves):
        print >> out, '  <taxon idref="taxon_%s"/>' % N[leaf]
    print >> out, '</taxa>'
    return out.getvalue().rstrip()

def get_mrca_stat_defn(name):
    """
    Get the second chunk of http://beast.bio.ed.ac.uk/Tutorial_3.1
    @param name: the name of the internal vertex
    """
    out = StringIO()
    print >> out, '<tmrcaStatistic id="tmrca_subset_%s" name="tmrca_%s">' % (
            name, name)
    print >> out, '  <treeModel idref="treeModel1"/>'
    print >> out, '  <mrca>'
    print >> out, '    <taxa idref="taxa_subset_%s"/>' % name
    print >> out, '  </mrca>'
    print >> out, '</tmrcaStatistic>'
    return out.getvalue().rstrip()

def get_starting_tree_defn(R, B, N_leaves):
    """
    """
    out = StringIO()
    N_aug = dict((v, 'taxon_' + name) for v, name in N_leaves.items())
    print >> out, '<newick id="startingTree">'
    print >> out, FtreeIO.RBN_to_newick(R, B, N_aug)
    print >> out, '</newick>'
    return out.getvalue().rstrip()

def get_leaf_taxon_defn(leaf_names):
    out = StringIO()
    print >> out, '<taxa id="taxa1">'
    for name in leaf_names:
        print >> out, '  <taxon id="taxon_%s">' % name
        print >> out, '    <date value="0" units="years" direction="forwards"/>'
        print >> out, '  </taxon>'
    print >> out, '</taxa>'
    return out.getvalue().rstrip()

def sample_b_to_rate(R, epsrate):
    v_to_source = R_to_v_to_source(R)
    for v in R_to_preorder(R):
        p = v_to_source.get(v, None)


def sample_v_to_seq(R, B, nsites):
    """
    @param R: directed topology
    @param B: branch lengths in expected substitutions
    @param nsites: sample this many independent sites of the ungapped alignment
    @return: map from vertex to sequence
    """
    v_to_seq = defaultdict(list)
    for i in range(nsites):
        for v, nt in sample_jc_column(R, B).items():
            v_to_seq[v].append(nt)
    return v_to_seq


def get_response_content(fs):
    # init the response and get the user variables
    out = StringIO()
    nleaves = fs.nleaves
    nsites = fs.nsites
    epsrate = fs.epsrate
    # sample the coalescent tree with timelike branch lengths
    R, B = kingman.sample(fs.nleaves)
    # get the leaf vertex names
    N = dict(zip(range(nleaves), string.uppercase[:nleaves]))
    N_leaves = dict(N)
    # get the internal vertex names
    v_to_leaves = R_to_v_to_leaves(R)
    for v, leaves in sorted(v_to_leaves.items()):
        if len(leaves) > 1:
            N[v] = ''.join(sorted(N[leaf] for leaf in leaves))
    # sample the alignment
    v_to_seq = {}
    #
    print >> out, '<?xml version="1.0"?>'
    print >> out, '<beast>'
    print >> out
    print >> out, '<!--'
    print >> out, 'predefine the taxa as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Introduction_to_XML_format'
    print >> out, '-->'
    print >> out, get_leaf_taxon_defn(list(string.uppercase[:nleaves]))
    print >> out
    print >> out, '<!--'
    print >> out, 'specify the starting tree as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Tutorial_4'
    print >> out, '-->'
    print >> out, get_starting_tree_defn(R, B, N_leaves)
    print >> out
    print >> out, '<!--'
    print >> out, 'connect the tree model as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Tutorial_4'
    print >> out, '-->'
    print >> out, g_tree_model_defn
    print >> out
    print >> out, '<!--'
    print >> out, 'create a list of taxa for which to constrain the mrca as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Tutorial_3.1'
    print >> out, '-->'
    for v, leaves in sorted(v_to_leaves.items()):
        if len(leaves) > 1:
            print >> out, get_mrca_subset_defn(N, v, leaves)
    print >> out
    print >> out, '<!--'
    print >> out, 'create a tmrcaStatistic that will record the height as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Tutorial_3.1'
    print >> out, '-->'
    for v, leaves in sorted(v_to_leaves.items()):
        if len(leaves) > 1:
            print >> out, get_mrca_stat_defn(N[v])
    print >> out
    print >> out, '</beast>'
    # return the response
    return out.getvalue()

