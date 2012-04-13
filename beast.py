"""
A module for interfacing with the BEAST MCMC java program.

This assumes that BEAST is on the path and that BEAGLE has been installed.
1) prepare() the base path
2) write the "<basepath>/myjob.xml" file referencing "myjob.log"
3) run_beast(basepath)
4) run_loganalyser(basepath)
"""

import subprocess
import time
import os
import tempfile

import Util


def prepare():
    """
    Create the job directory.
    @return: path to a newly created job directory
    """
    # reserve a name
    data = 'reserving a unique random local name...\n'
    tmp_path = Util.create_tmp_file(data=data, prefix='', suffix='')
    tmp_name = os.path.basename(tmp_path)
    t = time.gmtime()
    time_str = '%04d%02d%02d' % (t.tm_year, t.tm_mon, t.tm_mday)
    compound_str = '_'.join(('beast', time_str, tmp_name))
    new_dir = os.path.join(tempfile.gettempdir(), compound_str)
    # create the input and output directories
    os.makedirs(new_dir)
    # return the base directory
    return new_dir

def run_beast(base_path, jar_path):
    """
    The base path will be used as the working directory.
    Whatever miscellaneous intermediate files beast wants to make
    will be created in that directory.
    The base path is assumed to contain the xml named myjob.xml.
    The log file named myjob.log will be created in that directory.
    Also myjob-beast.out and myjob-beast.err.
    @param base_path: base path to run this job
    @param jar_path: path to the beast jar file
    @return: None
    """
    out_path = os.path.join(base_path, 'myjob-beast.out')
    err_path = os.path.join(base_path, 'myjob-beast.err')
    xml_path = os.path.join(base_path, 'myjob.xml')
    args = [
            #'beast',
            'java', '-jar', jar_path,
            '-beagle', '-beagle_CPU', '-beagle_SSE', '-beagle_double',
            '-working',
            xml_path]
    with open(out_path, 'w') as fout:
        with open(err_path, 'w') as ferr:
            p = subprocess.Popen(args, stdout=fout, stderr=ferr)
            p.communicate()

def run_loganalyser(base_path):
    """
    @param base_path: base path used in beast analysis
    @return: None
    """
    log_path = os.path.join(base_path, 'myjob.log')
    txt_path = os.path.join(base_path, 'myjob-loganalyser.txt')
    out_path = os.path.join(base_path, 'myjob-loganalyser.out')
    err_path = os.path.join(base_path, 'myjob-loganalyser.err')
    args = [
            'loganalyser',
            '-hpd', '-ess', log_path, txt_path]
    with open(out_path, 'w') as fout:
        with open(err_path, 'w') as ferr:
            p = subprocess.Popen(args, cwd=base_path, stdout=fout, stderr=ferr)
            p.communicate()

