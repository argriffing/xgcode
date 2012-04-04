"""
Sample some data so that BEAST can later guess the rate autocorrelation.

Autocorrelated rates on branches are sampled by assuming that
the implicit 'root branch' has rate 1.0
and that child branches have a rate that is sampled according
to an exponential distribution whose mean is the parent branch rate.
BEAST defines autocorrelation as the correlation between
parent lineage rates and child lineage rates.
In its covariance calculation, BEAST uses Bessel's correction
referred to in its comments as a correction so that the estimate is ML,
but this is a mistake because the correction is actually so that
the estimate is unbiased.
Also it has a function called covariance
that instead computes the correlation.
"""

from StringIO import StringIO
import string
import math
import random
from collections import defaultdict
import scipy
from scipy import optimize

import Form
import FormOut
import kingman
import FtreeIO
import Ftree
import RateMatrix
import MatrixUtil
import PhyLikelihood
import Fasta
import Newick

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

# mean param was originally  2.3E-4
# stdev param was originally 0.1
g_uncorrelated_relaxed_clock_info = """
<!--
the following block is lightly adapted from
http://beast.bio.ed.ac.uk/Hepatitis_C_subtype_1a,_genes_E1_and_E2
-->
<!-- The uncorrelated relaxed clock (Drummond, Ho, Phillips & Rambaut, 2006) -->
<discretizedBranchRates id="branchRates">
<treeModel idref="treeModel1"/>
<distribution>
<logNormalDistributionModel meanInRealSpace="true">
<mean>
  <parameter id="ucld.mean" value="1.0" lower="0.0" upper="100.0"/>
</mean>
<stdev>
  <parameter id="ucld.stdev" value="1.0" lower="0.0" upper="10.0"/>
</stdev>
</logNormalDistributionModel>
</distribution>
<rateCategories>
  <parameter id="branchRates.categories" dimension="120"/>
</rateCategories>
</discretizedBranchRates>
<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">
  <treeModel idref="treeModel1"/>
  <discretizedBranchRates idref="branchRates"/>
</rateStatistic>
<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">
  <treeModel idref="treeModel1"/>
  <discretizedBranchRates idref="branchRates"/>
</rateStatistic>
<rateCovarianceStatistic id="covariance" name="covariance">
  <treeModel idref="treeModel1"/>
  <discretizedBranchRates idref="branchRates"/>
</rateCovarianceStatistic>
""".strip()

g_likelihood_info = """
<patterns id="patterns1">
  <alignment idref="alignment1"/>
</patterns>
<!-- this is supposed to be Jukes-Cantor -->
<!-- see for example the question in the http://beast.bio.ed.ac.uk/FAQ -->
<gtrModel id="gtr1">
    <frequencies>
        <frequencyModel dataType="nucleotide">
        <!-- [preventing empirical freqs] alignment idref="alignment1"/ -->
        <frequencies>
            <!-- http://beast.bio.ed.ac.uk/Substitution_model_code -->
            <parameter id="gtr.frequencies" value="0.25 0.25 0.25 0.25"/>
        </frequencies>
        </frequencyModel>
    </frequencies>
    <rateAC> <parameter id="gtr.ac" value="1.0"/> </rateAC>
    <rateAG> <parameter id="gtr.ag" value="1.0"/> </rateAG>
    <rateAT> <parameter id="gtr.at" value="1.0"/> </rateAT>
    <rateCG> <parameter id="gtr.cg" value="1.0"/> </rateCG>
    <rateGT> <parameter id="gtr.gt" value="1.0"/> </rateGT>
</gtrModel>
<siteModel id="siteModel1">
  <substitutionModel>
    <gtrModel idref="gtr1"/>
  </substitutionModel>
</siteModel>
<treeLikelihood id="treeLikelihood">
  <patterns idref="patterns1"/>
  <treeModel idref="treeModel1"/>
  <siteModel idref="siteModel1"/>
  <discretizedBranchRates idref="branchRates"/>
</treeLikelihood>
""".strip()

g_log_info = """
<log id="fileLog" logEvery="10" fileName="experimental-beast.log">
    <posterior idref="posterior"/>
    <treeLikelihood idref="treeLikelihood"/>
    <parameter idref="ucld.mean"/>
    <parameter idref="ucld.stdev"/>
    <rateStatistic idref="coefficientOfVariation"/>
    <rateCovarianceStatistic idref="covariance"/>
    <parameter idref="treeModel1.rootHeight"/>
</log>
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
            #Form.Float('epsrate', 'larger values mean more rate variation',
                #0.02, low_exclusive=0, high_exclusive=1),
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
    v_to_source = Ftree.R_to_v_to_source(R)
    for v in Ftree.R_to_preorder(R):
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
        print >> out, '  <taxon id="taxon_%s"/>' % name
        #print >> out, '    <date value="0" units="years"'
        #print >> out, '      direction="forwards"/>'
        #print >> out, '  </taxon>'
    print >> out, '</taxa>'
    return out.getvalue().rstrip()

def sample_b_to_rate(R):
    """
    The b in this function name means branch.
    @param R: directed topology
    @return: a sampled map from vertex to expected rate
    """
    b_to_rate = {}
    v_to_source = Ftree.R_to_v_to_source(R)
    for v in Ftree.R_to_preorder(R):
        p = v_to_source.get(v, None)
        if p is None:
            continue
        # sample a coefficient regardless of whether we use it
        # this is an obsolete method
        #log_coeff = (random.random() - 0.5) * epsrate
        #coeff = math.exp(log_coeff)
        curr_branch = frozenset([v, p])
        gp = v_to_source.get(p, None)
        if gp is None:
            parent_rate = 1.0
        else:
            prev_branch = frozenset([p, gp])
            parent_rate = b_to_rate[prev_branch]
        b_to_rate[curr_branch] = random.expovariate(1/parent_rate)
    return b_to_rate

def get_correlation(R, b_to_rate):
    """
    This tries to exactly replicate the BEAST statistic.
    """
    X = []
    Y = []
    v_to_source = Ftree.R_to_v_to_source(R)
    for p, v in R:
        gp = v_to_source.get(p, None)
        if gp is not None:
            X.append(b_to_rate[frozenset([gp, p])])
            Y.append(b_to_rate[frozenset([p, v])])
    xbar = sum(X) / len(X)
    ybar = sum(Y) / len(Y)
    xvar = sum((x - xbar)**2 for x in X) / (len(X) - 1)
    yvar = sum((y - ybar)**2 for y in Y) / (len(Y) - 1)
    xstd = math.sqrt(xvar)
    ystd = math.sqrt(yvar)
    xycorr_num = sum((x - xbar) * (y - ybar) for x, y in zip(X, Y))
    xycorr_den = xstd * ystd * len(zip(X, Y))
    xycorr = xycorr_num / xycorr_den
    return xycorr


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

def get_alignment_defn(leaves, N, v_to_seq):
    """
    @param leaves: leaf vertices
    @param N: map from vertex to name
    @param v_to_seq: map from vertex to sequence
    """
    out = StringIO()
    print >> out, '<alignment id="alignment1" dataType="nucleotide">'
    for v in sorted(leaves):
        print >> out, '  <sequence>'
        print >> out, '    <taxon idref="taxon_%s"/>' % N[v]
        print >> out, '    ' + ''.join(v_to_seq[v])
        print >> out, '  </sequence>'
    print >> out, '</alignment>'
    return out.getvalue().rstrip()

def get_mcmc_defn(v_to_leaves, v_to_age, N):
    out = StringIO()
    print >> out, '<operators id="operators">'
    print >> out, '  <scaleOperator scaleFactor="0.75" weight="3">'
    print >> out, '    <parameter idref="ucld.stdev"/>'
    print >> out, '  </scaleOperator>'
    print >> out, '  <scaleOperator scaleFactor="0.75" weight="3">'
    print >> out, '    <parameter idref="ucld.mean"/>'
    print >> out, '  </scaleOperator>'
    #print >> out, '  <scaleOperator scaleFactor="0.75" weight="30"'
    #print >> out, '      autoOptimize="false">'
    #print >> out, '    <parameter idref="branchRates.categories"/>'
    #print >> out, '  </scaleOperator>'
    print >> out, '  <swapOperator size="1" weight="30" autoOptimize="false">'
    print >> out, '    <compoundParameter idref="branchRates.categories"/>'
    print >> out, '  </swapOperator>'
    print >> out, '</operators>'
    print >> out, '<mcmc id="mcmc1" chainLength="1000" autoOptimize="true">'
    print >> out, '  <posterior id="posterior">'
    """
    print >> out, '  <prior id="prior">'
    for v, leaves in sorted(v_to_leaves.items()):
        if len(leaves) != 1:
            print >> out, '    <uniformPrior lower="%s" upper="%s">' % (
                    v_to_age[v] * .99, v_to_age[v] * 1.01)
            print >> out, '    <tmrcaStatistic idref="tmrca_subset_%s"/>' % N[v]
            print >> out, '    </uniformPrior>'
    print >> out, '  </prior>'
    """
    print >> out, '  <likelihood id="likelihood">'
    print >> out, '    <treeLikelihood idref="treeLikelihood"/>'
    print >> out, '  </likelihood>'
    print >> out, '  </posterior>'
    print >> out, '  <operators idref="operators"/>'
    print >> out, g_log_info
    print >> out, '</mcmc>'
    return out.getvalue().rstrip()

class Opt:
    """
    This is for maximum likelihood search of lineage rates.
    """
    def __init__(self, R, B, N_leaves, alignment):
        """
        The vertices should be consecutive integers starting at zero.
        The largest vertex should be the root.
        """
        self.R = R
        self.B = B
        self.N_leaves = N_leaves
        self.alignment = alignment
    def __call__(self, X_logs):
        """
        The vth entry of X corresponds to the log rate of the branch above v.
        Return the quantity to be minimized (the neg log likelihood).
        @param X: vector of branch rate logs
        @return: negative log likelihood
        """
        X = [math.exp(x) for x in X_logs]
        B_subs = {}
        for v_parent, v_child in self.R:
            edge = frozenset([v_parent, v_child])
            r = X[v_child]
            t = self.B[edge]
            B_subs[edge] = r * t
        newick_string = FtreeIO.RBN_to_newick(self.R, B_subs, self.N_leaves)
        tree = Newick.parse(newick_string, Newick.NewickTree)
        # define the rate matrix object; horrible
        dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix() 
        ordered_states = list('ACGT') 
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(
                dictionary_rate_matrix, ordered_states, ordered_states)
        rate_matrix_object = RateMatrix.RateMatrix(
                row_major_rate_matrix, ordered_states) 
        # get the log likelihood
        ll = PhyLikelihood.get_log_likelihood(
                tree, self.alignment, rate_matrix_object)
        return -ll



def get_response_content(fs):
    # init the response and get the user variables
    out = StringIO()
    nleaves = fs.nleaves
    nvertices = nleaves * 2 - 1
    nbranches = nvertices - 1
    nsites = fs.nsites
    # sample the coalescent tree with timelike branch lengths
    R, B = kingman.sample(fs.nleaves)
    r = Ftree.R_to_root(R)
    # get the leaf vertex names
    N = dict(zip(range(nleaves), string.uppercase[:nleaves]))
    N_leaves = dict(N)
    # get the internal vertex names
    v_to_leaves = R_to_v_to_leaves(R)
    for v, leaves in sorted(v_to_leaves.items()):
        if len(leaves) > 1:
            N[v] = ''.join(sorted(N[leaf] for leaf in leaves))
    # get vertex ages
    v_to_age = kingman.RB_to_v_to_age(R, B)
    # sample the rates on the branches
    b_to_rate = sample_b_to_rate(R)
    xycorr = get_correlation(R, b_to_rate)
    # define B_subs in terms of substitutions instead of time
    B_subs = dict((p, t * b_to_rate[p]) for p, t in B.items())
    # sample the alignment
    v_to_seq = sample_v_to_seq(R, B_subs, nsites)
    # get the log likelihood; this is kind of horrible
    pairs = [(N[v], ''.join(v_to_seq[v])) for v in range(nleaves)]
    headers, sequences = zip(*pairs)
    alignment = Fasta.create_alignment(headers, sequences)
    newick_string = FtreeIO.RBN_to_newick(R, B_subs, N_leaves)
    tree = Newick.parse(newick_string, Newick.NewickTree)
    dictionary_rate_matrix = RateMatrix.get_jukes_cantor_rate_matrix() 
    ordered_states = list('ACGT') 
    row_major_rate_matrix = MatrixUtil.dict_to_row_major(
            dictionary_rate_matrix, ordered_states, ordered_states)
    rate_matrix_object = RateMatrix.RateMatrix(
            row_major_rate_matrix, ordered_states) 
    ll = PhyLikelihood.get_log_likelihood(
            tree, alignment, rate_matrix_object)
    # get ll when rates are all 1.0
    newick_string = FtreeIO.RBN_to_newick(R, B, N_leaves)
    tree = Newick.parse(newick_string, Newick.NewickTree)
    ll_unity = PhyLikelihood.get_log_likelihood(
            tree, alignment, rate_matrix_object)
    # get ll when rates are numerically optimized
    # TODO incorporate the result into the xml file
    # TODO speed up the likelihood evaluation (beagle? C module?)
    #f = Opt(R, B, N_leaves, alignment)
    #X_logs = [0.0] * nbranches
    #result = scipy.optimize.fmin(f, X_logs, full_output=True)
    #print result
    #
    print >> out, '<?xml version="1.0"?>'
    print >> out, '<beast>'
    print >> out
    print >> out, '<!-- actual rate autocorrelation', xycorr, '-->'
    print >> out, '<!-- actual root height', v_to_age[r], '-->'
    print >> out, '<!-- actual log likelihood', ll, '-->'
    print >> out, '<!-- ll if rates were unity', ll_unity, '-->'
    print >> out
    print >> out, '<!--'
    print >> out, 'predefine the taxa as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Introduction_to_XML_format'
    print >> out, '-->'
    print >> out, get_leaf_taxon_defn(list(string.uppercase[:nleaves]))
    print >> out
    print >> out, '<!--'
    print >> out, 'define the alignment as in'
    print >> out, 'http://beast.bio.ed.ac.uk/Introduction_to_XML_format'
    print >> out, '-->'
    print >> out, get_alignment_defn(leaves, N, v_to_seq)
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
    print >> out, g_uncorrelated_relaxed_clock_info
    print >> out
    """
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
    """
    print >> out
    print >> out, g_likelihood_info
    print >> out
    print >> out, '<!--'
    print >> out, 'run the mcmc'
    print >> out, 'http://beast.bio.ed.ac.uk/Tutorial_3.1'
    print >> out, '-->'
    print >> out, get_mcmc_defn(v_to_leaves, v_to_age, N)
    print >> out
    print >> out, '</beast>'
    # return the response
    return out.getvalue()

