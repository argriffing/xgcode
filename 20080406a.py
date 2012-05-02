"""
Use PAML to estimate HKY85 parameters given nexus data.

The nexus data should have a tree and an alignment.
"""

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
import Util

# TODO use subprocess instead of popen2 to run paml
# and use the subprocess option to run in different working directory
# This whole module is kind of ancient and bad.

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

def run_in_tmp(fs, nexus):
    """
    @param fs: fieldstorage-like user options
    @param nexus: the nexus object that defines the tree and alignment
    @return: from_paml
    """
    # create the control object
    config = Paml.PamlConfig()
    config.set_hky()
    # create the baseml.ctl control file
    ctl_string = config.to_ctl_string()
    with open(Paml.baseml_ctl, 'wt') as fout:
        print >> fout, ctl_string
    # create the nexus object that defines the tree and alignment
    nexus = Nexus.get_sample_nexus_object()
    # create the baseml.newick tree file
    with open(Paml.baseml_newick, 'wt') as fout:
        print >> fout, nexus.tree.get_newick_string()
    # create the baseml.phylip alignment file
    s_phylip = Phylip.get_alignment_string_non_interleaved(nexus.alignment)
    with open(Paml.baseml_phylip, 'wt') as fout:
        print >> fout, s_phylip
    # run PAML
    args = [Config.baseml_exe_path, Paml.baseml_ctl]
    from_paml, to_paml = popen2.popen4(args)
    return from_paml

def get_response_content(fs):
    # read the nexus data
    nexus = Nexus.Nexus()
    try:
        nexus.load(StringIO(fs.nexus))
    except Nexus.NexusError as e:
        raise HandlingError(e)
    # run paml in a temp directory
    with Util.remember_cwd():
        os.chdir(Config.data_path)
        from_paml = run_in_tmp(fs, nexus)
    # get the paml output
    paml_debug_output = from_paml.read()
    # get the paml output file
    with open(Paml.baseml_out) as fin:
        d = Paml.parse_hky_output(fin)
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
