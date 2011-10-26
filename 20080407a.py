"""Use HyPhy to estimate HKY85 parameters given nexus data.

The nexus data should have a tree and an alignment.
"""

import math
from StringIO import StringIO
import os
import subprocess

from SnippetUtil import HandlingError
import Config
import Hyphy
import Nexus
import Form
import FormOut

# define the input nexus path
hyphy_nexus = os.path.join(Config.data_path, 'hyphy.nex')
# define the input batch file model path
hyphy_bf = os.path.join(Config.data_path, 'model.bf')
# tell hyphy to use this model
hky_hyphy_model_sio = StringIO()
print >> hky_hyphy_model_sio, 'VERBOSITY_LEVEL = 1;'
print >> hky_hyphy_model_sio, 'ACCEPT_BRANCH_LENGTHS = 1;'
print >> hky_hyphy_model_sio, 'DataSet spectrinData = ReadDataFile ("%s");' % hyphy_nexus
print >> hky_hyphy_model_sio, r"""
DataSetFilter filteredData = CreateFilter (spectrinData,1);
HarvestFrequencies (observedFreqs, filteredData, 1, 1, 1);
/* stationary distribution in nucleotide frequencies */
global  eqFreqA = 0.25;
global  eqFreqC = 0.25;
global  eqFreqG = 0.25;
global  eqFreqT := 1.0 - eqFreqA - eqFreqC - eqFreqG;   eqFreqT :> 0;
/* transition/transversion ratio */
global  kappa = 2.0;    kappa :> 0;
HKY85RateMatrix =
  {{*,a,kappa*a,a}
   {a,*,a,kappa*a}
   {kappa*a,a,*,a}
   {a,kappa*a,a,*}};
/* construct vectors from global variables */
estFreqs =
  {{eqFreqA,
    eqFreqC,
    eqFreqG,
    eqFreqT}
  };
/* assign models to trees */
Model m1 = (HKY85RateMatrix, estFreqs);
Tree givenTree = DATAFILE_TREE;
branchNames = BranchName (givenTree, -1);  /* vector populated with names for each branch in tree */
branchLengths = BranchLength (givenTree, -1);  /* vector populated with branch lengths */
for (k = 0; k < Columns(branchNames)-1; k = k+1)
{ 
  ExecuteCommands("givenTree."+branchNames[k]+".a:="+branchLengths[k]+";");
}
LikelihoodFunction lf = (filteredData, givenTree, "Log(SITE_LIKELIHOOD[0])");
Optimize(res, lf);
fprintf(stdout, "\n");
fprintf(stdout, res);
fprintf(stdout, lf);
""".strip()
hky_hyphy_model = hky_hyphy_model_sio.getvalue().strip()

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('nexus',
                'nexus data', Nexus.nexus_sample_string.strip()),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('outdebug', 'show debug info'),
                Form.CheckItem('outmodel', 'show the model')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the nexus data
    nexus = Nexus.Nexus()
    try:
        nexus.load(StringIO(fs.nexus))
    except Nexus.NexusError, e:
        raise HandlingError(e)
    # move to the data directory
    original_directory = os.getcwd()
    os.chdir(Config.data_path)
    # create the batch file
    fout = open(hyphy_bf, 'wt')
    print >> fout, hky_hyphy_model
    fout.close()
    # create the nexus file
    fout = open(hyphy_nexus, 'wt')
    print >> fout, nexus
    fout.close()
    # run hyphy
    cmd = [Config.hyphy_exe_path, hyphy_bf]
    p = subprocess.Popen(cmd, close_fds=True, stdout=subprocess.PIPE)
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
    print >> out, 'tree :', ns.givenTree.get_newick_string()
    print >> out, ''
    print >> out, 'log likelihood :', ns.lnL
    print >> out, ''
    print >> out, 'kappa :', ns.kappa
    print >> out, ''
    print >> out, 'A :', ns.eqFreqA
    print >> out, 'C :', ns.eqFreqC
    print >> out, 'G :', ns.eqFreqG
    print >> out, 'T :', ns.eqFreqT
    # return the results
    return out.getvalue()

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
