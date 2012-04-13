"""
Use BEAST to analyze an interval of the primate tutorial alignment.

Web service output is the loganalyser output.
The extended command line functionality
will eventually include
the analysis of multiple intervals
using a separate BEAST instance for each interval.
"""

from StringIO import StringIO
import argparse
import multiprocessing
import os

import Form
import FormOut
import beasttut
import beast

g_beast_jar_path = os.path.expanduser('~/BEASTv1.7.1/lib/beast.jar')

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('start',
                'sub-sequence start position (1-%d)' % beasttut.g_nchar,
                1, low=1, high=beasttut.g_nchar),
            Form.Integer('stop',
                'sub-sequence stop position (1-%d)' % beasttut.g_nchar,
                beasttut.g_nchar, low=1, high=beasttut.g_nchar),
            Form.Integer('nsamples', 'mcmc chain steps',
                8000, low=80, high=8000)]
    return form_objects

def get_form_out():
    return FormOut.Report('summary.txt')

def get_response_content(fs):
    # init the response and get the user variables
    start_pos = fs.start
    stop_pos = fs.stop
    nsamples = fs.nsamples
    out = StringIO()
    # get the xml contents
    xml_string = beasttut.get_xml_string(
            start_pos, stop_pos, nsamples, 'myjob.log',
            beasttut.get_header_seq_pairs())
    # prepare the base path for the beast analysis
    basepath = beast.prepare()
    with open(os.path.join(basepath, 'myjob.xml'), 'w') as fout:
        fout.write(xml_string + '\n')
    beast.run_beast(basepath, g_beast_jar_path)
    beast.run_loganalyser(basepath)
    # return the analysis
    with open(os.path.join(basepath, 'myjob-loganalyser.txt')) as fin:
        analysis_text = fin.read()
    # return the response
    return analysis_text

