"""Use HyPhy to estimate kappa and a frequency mixture given nexus data.

The nexus data should have a tree and an alignment.
"""

import math
from StringIO import StringIO
import os
import subprocess

from SnippetUtil import HandlingError
import Config
import Hyphy
import Newick
import RateMatrix
import SubModel
import PhyLikelihood
import Nexus
import Form
import FormOut

# define the input nexus path
hyphy_nexus = os.path.join(Config.data_path, 'hyphy.nex')
# define the input batch file model path
hyphy_bf = os.path.join(Config.data_path, 'model.bf')

def get_form():
    """
    @return: the body of a form
    """
    # define the default nexus string
    tree = get_sample_tree()
    mixture_model = get_sample_mixture_model()
    ncols = 200
    seed = 314159
    alignment = PhyLikelihood.simulate_alignment(
            tree, mixture_model, ncols, seed)
    nexus = Nexus.Nexus()
    nexus.tree = tree
    nexus.alignment = alignment
    nexus_string = str(nexus)
    # define the form objects
    form_objects = [
            Form.MultiLine('nexus', 'nexus data', nexus_string),
            Form.Integer('ncategories', 'use this many categories',
                3, low=1, high=5),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('outdebug', 'show debug info'),
                Form.CheckItem('outmodel', 'show the model'),
                Form.CheckItem('outcheck', 'show the likelihood and rates',
                    True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the nexus data
    nexus = Nexus.Nexus()
    try:
        nexus.load(StringIO(fs.nexus))
    except Nexus.NexusError as e:
        raise HandlingError(e)
    # move to the data directory
    original_directory = os.getcwd()
    os.chdir(Config.data_path)
    # create the batch file
    fout = open(hyphy_bf, 'wt')
    category_suffixes = [str(category+1) for category in range(fs.ncategories)]
    hky_hyphy_model = get_hyphy_model_string(hyphy_nexus, fs.ncategories)
    print >> fout, hky_hyphy_model 
    fout.close()
    # create the nexus file
    fout = open(hyphy_nexus, 'wt')
    print >> fout, nexus
    fout.close()
    # run hyphy
    p = subprocess.Popen([Config.hyphy_exe_path, hyphy_bf],
            close_fds=True, stdout=subprocess.PIPE)
    hyphy_output = p.stdout.read()
    # move back to the original directory
    os.chdir(original_directory)
    # read the hyphy output
    ns = Hyphy.get_hyphy_namespace(StringIO(hyphy_output))
    out = StringIO()
    if fs.outdebug:
        print >> out, get_hyphy_debug_info(hyphy_output)
        print >> out, ''
        print >> out, ''
    if fs.outmodel:
        print >> out, 'hyphy model:'
        print >> out, '---------------------------------------'
        print >> out, hky_hyphy_model
        print >> out, '---------------------------------------'
        print >> out, ''
        print >> out, ''
    if True:
        print >> out, 'reformatted hyphy output:'
        print >> out, '---------------------------------------'
        # show the log likelihood
        print >> out, 'log likelihood :', ns.lnL
        print >> out, ''
        # show the kappa value
        print >> out, 'kappa :', ns.kappa
        print >> out, ''
        category_blocks = []
        for suffix in category_suffixes:
            block = StringIO()
            print >> block, 'mixing proportion :', getattr(ns, 'catFreq'+suffix)
            print >> block, 'tree :', getattr(ns, 'tree'+suffix).get_newick_string()
            for nt in list('ACGT'):
                print >> block, nt, ':', getattr(ns, 'eqFreq'+nt+suffix)
            category_blocks.append(block.getvalue().strip())
        print >> out, '\n\n'.join(category_blocks)
        print >> out, '---------------------------------------'
        print >> out, ''
        print >> out, ''
    if fs.outcheck:
        # get the raw matrices
        matrices = []
        for suffix in category_suffixes:
            nt_dict = {}
            for nt in list('ACGT'):
                nt_dict[nt] = getattr(ns, 'eqFreq'+nt+suffix)
            total = float(sum(nt_dict.values()))
            nt_dict = dict((k, v/total) for k, v in nt_dict.items())
            matrix = RateMatrix.get_unscaled_hky85_rate_matrix(
                    nt_dict, ns.kappa)
            matrices.append(matrix)
        raw_matrix_rates = [matrix.get_expected_rate() for matrix in matrices]
        category_weights = []
        for suffix in category_suffixes:
            category_weights.append(getattr(ns, 'catFreq'+suffix))
        total = float(sum(category_weights))
        category_distribution = [weight / total for weight in category_weights]
        mixture_model = SubModel.MixtureModel(category_distribution, matrices)
        raw_mixture_rate = mixture_model.get_expected_rate()
        # rescale the mixture model
        # 0.75 is the expected rate of the initial model
        r1 = 0.75
        scaling_factor = r1
        mixture_model.rescale(scaling_factor)
        recomputed_log_likelihood = PhyLikelihood.get_log_likelihood(
                nexus.tree, nexus.alignment, mixture_model)
        print >> out, 'recomputed likelihood and rates:'
        print >> out, '---------------------------------------'
        print >> out, 'log likelihood :', recomputed_log_likelihood
        print >> out, ''
        print >> out, 'rate :', raw_mixture_rate
        print >> out, ''
        for rate, suffix in zip(raw_matrix_rates, category_suffixes):
            print >> out, 'rate%s : %s' % (suffix, rate)
        print >> out, '---------------------------------------'
        print >> out, ''
        print >> out, ''
    # return the response
    return out.getvalue()

def get_hyphy_model_string(nexus_path, ncategories):
    """
    @param nexus_path: the path to a nexus file with a tree and alignment
    @param ncategories: the number of categories in the mixture
    @return: a hyphy batch file
    """
    out = StringIO()
    print >> out, 'VERBOSITY_LEVEL = 1;'
    print >> out, 'ACCEPT_BRANCH_LENGTHS = 1;'
    print >> out, 'DataSet spectrinData = ReadDataFile ("%s");' % nexus_path
    print >> out, 'DataSetFilter filteredData = CreateFilter (spectrinData,1);'
    print >> out, 'HarvestFrequencies (observedFreqs, filteredData, 1, 1, 1);'
    print >> out, '/* stationary distribution of nucleotide frequencies */'
    for category in range(ncategories):
        suffix = str(category + 1)
        print >> out, 'global  eqFreqA%s = 0.25;' % suffix
        print >> out, 'global  eqFreqC%s = 0.25;' % suffix
        print >> out, 'global  eqFreqG%s = 0.25;' % suffix
        print >> out, 'global  eqFreqT%s := 1.0 - eqFreqA%s - eqFreqC%s - eqFreqG%s;' % tuple([suffix]*4)
        print >> out, 'eqFreqT%s :> 0;' % suffix
    print >> out, '/* transition/transversion ratio */'
    print >> out, 'global  kappa = 1.0;'
    print >> out, 'kappa :> 0;'
    print >> out, '/* hky rate matrix */'
    print >> out, 'HKY85RateMatrix ='
    print >> out, '  {{*,a,kappa*a,a}'
    print >> out, '   {a,*,a,kappa*a}'
    print >> out, '   {kappa*a,a,*,a}'
    print >> out, '   {a,kappa*a,a,*}'
    print >> out, '  };'
    print >> out, '/* construct vectors from global variables */'
    for category in range(ncategories):
        suffix = str(category + 1)
        print >> out, 'estFreqs%s =' % suffix
        print >> out, '  {{eqFreqA%s,' % suffix
        print >> out, '    eqFreqC%s,' % suffix
        print >> out, '    eqFreqG%s,' % suffix
        print >> out, '    eqFreqT%s}' % suffix
        print >> out, '  };'
    print >> out, '/* assign models to trees */'
    for category in range(ncategories):
        suffix = str(category + 1)
        print >> out, 'Model m%s = (HKY85RateMatrix, estFreqs%s);' % (suffix, suffix)
        print >> out, 'Tree tree%s = DATAFILE_TREE;' % suffix
        print >> out, '/* vector populated with names for each branch in tree */'
        print >> out, 'branchNames = BranchName (tree%s, -1);' % suffix
        print >> out, '/* vector populated with branch lengths */'
        print >> out, 'branchLengths = BranchLength (tree%s, -1);' % suffix
        print >> out, 'for (k = 0; k < Columns(branchNames)-1; k = k+1)'
        print >> out, '{ '
        print >> out, '  ExecuteCommands("tree%s."+branchNames[k]+".a:="+branchLengths[k]+";");' % suffix
        print >> out, '}'
    print >> out, '/* define the mixing proportions and give them initial values */'
    for category in range(ncategories - 1):
        suffix = str(category + 1)
        print >> out, 'global catFreq%s = %f;' % (suffix, 1.0 / ncategories)
    cat_freqs = ['catFreq%d' % (i+1) for i in range(ncategories)]
    rhs_string = ' - '.join(['1.0'] + cat_freqs[:-1])
    print >> out, 'global catFreq%d := %s;' % (ncategories, rhs_string)
    print >> out, 'catFreq%d :> 0;' % (ncategories)
    print >> out, '/* construct the likelihood function */'
    likelihood_function_arguments = []
    for category in range(ncategories):
        suffix = str(category + 1)
        likelihood_function_arguments.append('filteredData')
        likelihood_function_arguments.append('tree' + suffix)
    likelihood_terms = []
    for category in range(ncategories):
        likelihood_terms.append('SITE_LIKELIHOOD[%d]*catFreq%d' % (category, category+1))
    likelihood_function_arguments.append('"Log(%s)"' % ' + '.join(likelihood_terms))
    print >> out, 'LikelihoodFunction lf = (%s);' % ', '.join(likelihood_function_arguments)
    print >> out, '/* search for maximum likelihood parameter estimates */'
    print >> out, 'Optimize(res, lf);'
    print >> out, '/* show the results */'
    print >> out, r'fprintf(stdout, "\n");'
    print >> out, 'fprintf(stdout, res);'
    print >> out, 'fprintf(stdout, lf);'
    return out.getvalue().strip()

def get_hyphy_debug_info(hyphy_output):
    """
    @param hyphy_output: the string representing the hyphy output
    @return: a string explaining how the output was interpreted
    """
    ns = Hyphy.get_hyphy_namespace(StringIO(hyphy_output))
    out = StringIO()
    print >> out, 'raw hyphy output:'
    print >> out, '---------------------------------------'
    print >> out, hyphy_output
    print >> out, '---------------------------------------'
    print >> out, ''
    print >> out, ''
    print >> out, 'processed hyphy output lines:'
    print >> out, '---------------------------------------'
    for i, line in enumerate(ns.get_processed_lines()):
        print >> out, i, ':', line
    print >> out, '---------------------------------------'
    print >> out, ''
    print >> out, ''
    print >> out, 'hyphy namespace object dictionary:'
    print >> out, '---------------------------------------'
    print >> out, ns.__dict__
    print >> out, '---------------------------------------'
    return out.getvalue().strip()

def get_sample_mixture_model():
    """
    @return: a mixture model that is used to generate the default nexus data
    """
    # define the model
    kappa = 2
    category_distribution = [.1, .4, .5]
    nt_dicts = [
            {'A' : .1, 'C' : .4, 'G' : .4, 'T' : .1},
            {'A' : .2, 'C' : .3, 'G' : .3, 'T' : .2},
            {'A' : .25, 'C' : .25, 'G' : .25, 'T' : .25}
            ]
    # create a mixture model from the variables that define the model
    rate_matrix_objects = []
    for nt_dict in nt_dicts:
        rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(
                nt_dict, kappa)
        rate_matrix_objects.append(rate_matrix_object)
    mixture_model = SubModel.MixtureModel(
            category_distribution, rate_matrix_objects)
    mixture_model.normalize()
    return mixture_model

def get_sample_tree():
    """
    @return: a mixture model that is used to generate the default nexus data
    """
    tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    return tree
