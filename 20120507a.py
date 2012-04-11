"""
Run BEAST on several sub-alignments. [UNFINISHED]

Plot several HPD intervals for a few statistics.
Colors are according to interval width (sub-alignment length).
"""

from StringIO import StringIO
import os
import subprocess
import argparse
import logging
import multiprocessing

import Form
import FormOut
import const
import Util
import RUtil
import hpcutil
import beasttut


g_nchar = 898
g_beast_root = os.path.expanduser('~/svn-repos/beast-mcmc-read-only')

g_start_stop_pairs = (
        # 8 of length 57
        (1, 57),
        (58, 114),
        (1 + 57*2, 57*3),
        (1 + 57*3, 57*4),
        (1 + 57*4, 57*5),
        (1 + 57*5, 57*6),
        (1 + 57*6, 57*7),
        (400, 456),
        # 4 of length 57*2
        (1, 114),
        (115, 228),
        (1 + 57*2*2, 57*2*3),
        (1 + 57*2*3, 57*2*4),
        # 2 of length 57*2*2
        (1, 228),
        (229, 456),
        # 1 of length 57*2*2*2
        (1, 456),
        )

g_headers = (
        'sequence.length',
        'midpoint',
        'mean.low', 'mean.mean', 'mean.high',
        'var.low', 'var.mean', 'var.high',
        'cov.low' ,'cov.mean', 'cov.high')


class RemoteBeast(hpcutil.RemoteBrc):
    def __init__(self, start_stop_pairs, nsamples):
        hpcutil.RemoteBrc.__init__(self)
        self.start_stop_pairs = start_stop_pairs
        self.nsamples = nsamples
        self.remote_beast_sh_path = os.path.join(
                '/brc_share/brc/argriffi/packages/BEASTv1.7.1',
                'bin/beast')
        self.remote_beast_jar_path = os.path.join(
                '/brc_share/brc/argriffi/packages/BEASTv1.7.1',
                'lib/beast.jar')
        self.local_log_paths = None
    def preprocess(self):
        # for each (start_pos, stop_pos) pair
        # define a log filename and
        # create a beast xml and a bsub script
        self.local_log_paths = []
        for i, start_stop_pair in enumerate(self.start_stop_pairs):
            start_pos, stop_pos = start_stop_pair
            # define the log file paths
            log_name = 'primate-tut-%d-%d.log' % start_stop_pair
            local_log_path = os.path.join(
                    self.local.get_out(), log_name)
            remote_log_path = os.path.join(
                    self.remote.get_out(), log_name)
            self.local_log_paths.append(local_log_path)
            # define the xml file paths
            xml_name = 'primate-tut-%d-%d.xml' % start_stop_pair
            local_xml_path = os.path.join(
                    self.local.get_in_contents(), xml_name)
            remote_xml_path = os.path.join(
                    self.remote.get_in_contents(), xml_name)
            # define the local bsub path
            bsub_name = 'primate-tut-%d-%d.bsub' % start_stop_pair
            local_bsub_path = os.path.join(
                    self.local.get_in_bsubs(), bsub_name)
            # create the xml file
            with open(local_xml_path, 'w') as fout:
                print >> fout, beasttut.get_xml_string(
                        start_pos, stop_pos, self.nsamples, remote_log_path)
            # create the bsub file
            with open(local_bsub_path, 'w') as fout:
                stdout_path = os.path.join(self.remote.get_out(),
                        'out.tut-%d-%d' % start_stop_pair)
                stderr_path = os.path.join(self.remote.get_out(),
                        'err.tut-%d-%d' % start_stop_pair)
                # name the job
                print >> fout, '#BSUB -J tut-%d-%d' % start_stop_pair
                # suggest the brc queue
                print >> fout, '#BSUB -q brc'
                # redirect stdout
                print >> fout, '#BSUB -o', stdout_path
                # redirect stderr
                print >> fout, '#BSUB -e', stderr_path
                # run BEAST using a java path suggested by Gary Howell
                print >> fout, '/usr/local/apps/java/jre1.6.0_31/bin/java',
                print >> fout, '-jar',
                print >> fout, self.remote_beast_jar_path, remote_xml_path

def make_xml(start_pos, stop_pos, nsamples):
    """
    This is for non-hpc only.
    @return: location of xml file, location of log file
    """
    log_loc = Util.get_tmp_filename(prefix='beast', suffix='.log')
    xml_string = beasttut.get_xml_string(
            start_pos, stop_pos, nsamples, log_loc)
    xml_loc = Util.create_tmp_file(xml_string, prefix='beast', suffix='.xml')
    return xml_loc, log_loc

def run_beast(xml_loc):
    """
    This is for non-hpc only.
    """
    args = (
            'java',
            '-jar',
            os.path.join(g_beast_root, 'build', 'dist', 'beast.jar'),
            xml_loc,
            )
    subprocess.call(args)
    

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('start',
                'sub-sequence start position (1-%d)' % g_nchar,
                1, low=1, high=g_nchar),
            Form.Integer('stop',
                'sub-sequence stop position (1-%d)' % g_nchar,
                g_nchar, low=1, high=g_nchar)]
    return form_objects

def get_form_out():
    return FormOut.Html()

def get_value_lists(start_pos, stop_pos, nsamples):
    """
    Command-line and also web based but not hpc-based.
    """
    # input validation
    if stop_pos < start_pos:
        raise ValueError('the stop pos must be after the start pos')
    # create the xml describing the analysis
    xml_loc, log_loc = make_xml(start_pos, stop_pos, nsamples)
    print 'log file location:', log_loc
    # run beast
    run_beast(xml_loc)
    # read the log file
    return beasttut.read_log(log_loc, nsamples)

def get_response_content(fs):
    # init the response and get the user variables
    start_pos = fs.start
    stop_pos = fs.stop
    nsamples = 8000
    out = StringIO()
    # do the analysis
    means, variations, covariances = get_value_lists(
            start_pos, stop_pos, nsamples)
    values_names_pairs = (
            (means, 'mean rate among branches'),
            (variations, 'coefficient of variation of rates among branches'),
            (covariances, 'correlation of parent and child branch rates'))
    print >> out, beasttut.get_html(values_names_pairs)
    # return the response
    return out.getvalue()

def main(args):
    # set up the logger
    f = logging.getLogger('toplevel.logger')
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter('%(message)s %(asctime)s'))
    f.addHandler(h)
    if args.verbose:
        f.setLevel(logging.DEBUG)
    else:
        f.setLevel(logging.WARNING)
    if args.run == 'hpc':
        r = RemoteBeast(g_start_stop_pairs, args.nsamples)
        f.info('run BEAST remotely')
        r.run(verbose=args.verbose)
        f.info('(local) build the R table string and scripts')
        table_string, scripts = beasttut.get_table_string_and_scripts_from_logs(
                g_start_stop_pairs, r.local_log_paths, args.nsamples)

    elif args.run == 'serial':
        f.info('(local) run BEAST serially locally and build the R stuff')
        table_string, scripts = beasttut.get_table_string_and_scripts(
                args.nsamples)
    elif args.run == 'parallel':
        f.info('(local) run BEAST locally in parallel and build the R stuff')
        table_string, scripts = beasttut.get_table_string_and_scripts_par(
                args.nsamples)
    else:
        raise ValueError('invalid execution model')
    f.info('(local) create the composite R script')
    out = StringIO()
    print >> out, 'library(ggplot2)'
    print >> out, 'par(mfrow=c(3,1))'
    for script in scripts:
        print >> out, script
    comboscript = out.getvalue()
    f.info('(local) run R to create the pdf')
    device_name = Form.g_imageformat_to_r_function['pdf']
    retcode, r_out, r_err, image_data = RUtil.run_plotter( 
        table_string, comboscript, device_name, keep_intermediate=True) 
    if retcode: 
        raise RUtil.RError(r_err) 
    f.info('(local) write the .pdf file')
    with open(args.outfile, 'wb') as fout:
        fout.write(image_data)
    f.info('(local) return from toplevel')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--outfile',
            default='beast-analysis.pdf',
            help='write this pdf file')
    parser.add_argument('--nsamples',
            default=8000, type=int,
            help='let the BEAST MCMC generate this many samples')
    parser.add_argument('--run', choices=('serial', 'parallel', 'hpc'),
            default='serial',
            help='execution model')
    parser.add_argument('-v', '--verbose',
            action='store_true',
            help='show more info')
    main(parser.parse_args())

