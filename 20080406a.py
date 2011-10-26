"""Use PAML to estimate HKY85 parameters given nexus data.

The nexus data should have a tree and an alignment.
"""

import math
from StringIO import StringIO
import os
import popen2

from SnippetUtil import HandlingError
import Nexus
import Paml
import Config
import Phylip
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('nexus', 'nexus data',
                Nexus.nexus_sample_string.strip()),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('outdebug', 'show debug info'),
                Form.CheckItem('outmodel', 'show the model')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # read the nexus data
    nexus = Nexus.Nexus()
    try:
        nexus.load(StringIO(fs.nexus))
    except Nexus.NexusError as e:
        raise HandlingError(e)
    # define some paths
    baseml_ctl = os.path.join(Config.data_path, 'baseml.ctl')
    baseml_phylip = os.path.join(Config.data_path, 'baseml.phylip')
    baseml_newick = os.path.join(Config.data_path, 'baseml.newick')
    baseml_out = os.path.join(Config.data_path, 'baseml.out')
    # change to the target directory
    original_directory = os.getcwd()
    os.chdir(Config.data_path)
    # create the baseml.ctl control file
    config = Paml.PamlConfig()
    config.set_hky()
    fout = open(baseml_ctl, 'wt')
    ctl_string = config.to_ctl_string()
    print >> fout, ctl_string
    fout.close()
    # create the nexus object that defines the tree and alignment
    nexus = Nexus.get_sample_nexus_object()
    # create the baseml.newick tree file
    fout = open(baseml_newick, 'wt')
    print >> fout, nexus.tree.get_newick_string()
    fout.close()
    # create the baseml.phylip alignment file
    fout = open(baseml_phylip, 'wt')
    s_phylip = Phylip.get_alignment_string_non_interleaved(nexus.alignment)
    print >> fout, s_phylip
    fout.close()
    # run PAML
    ctl_path = baseml_ctl
    from_paml, to_paml = popen2.popen4([Config.baseml_exe_path, ctl_path])
    # change back to the original directory
    os.chdir(original_directory)
    # get the paml output
    paml_debug_output = from_paml.read()
    # get the paml output file
    fin = open(Paml.baseml_out)
    d = Paml.parse_hky_output(fin)
    fin.close()
    out = StringIO()
    if fs.outdebug:
        print >> out, 'baseml stdout:'
        print >> out, '-----------------------------------------'
        print >> out, paml_debug_output.strip()
        print >> out, '-----------------------------------------'
        print >> out, ''
        print >> out, ''
    if fs.outmodel:
        print >> out, 'baseml control file:'
        print >> out, '-----------------------------------------'
        print >> out, ctl_string.strip()
        print >> out, '-----------------------------------------'
        print >> out, ''
        print >> out, ''
    print >> out, 'log likelihood :', d['lnL']
    print >> out, ''
    print >> out, 'kappa :', d['kappa']
    print >> out, ''
    for nt in list('ACGT'):
        print >> out, ('%s : ' % nt) + str(d[nt])
    # return the results
    return out.getvalue()
