"""
Use BEAST to analyze a column permutation of the primate tutorial alignment.

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
import logging

import Form
import FormOut
import beasttut
import beast
import RUtil

g_beast_jar_path = os.path.expanduser('~/BEASTv1.7.1/lib/beast.jar')

g_ncols_max = 456

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


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.Integer('ncolumns',
                'use this many columns (1-%d)' % g_ncols_max,
                g_ncols_max, low=1, high=g_ncols_max),
            Form.Integer('nsamples', 'mcmc chain steps',
                8000, low=80, high=8000)]
    return form_objects

def get_form_out():
    return FormOut.Report('summary')

def get_table_string_and_scripts(start_stop_pairs, nsamples, header_seq_pairs):
    """
    Command-line only.
    """
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_pos, stop_pos in start_stop_pairs:
        sequence_length = stop_pos - start_pos + 1
        midpoint = (start_pos + stop_pos) / 2.0
        arr = get_loganalysis_array(
                start_pos, stop_pos, nsamples, header_seq_pairs)
        mean_low = arr[1][4]
        mean_mean = arr[1][1]
        mean_high = arr[1][5]
        var_low = arr[2][4]
        var_mean = arr[2][1]
        var_high = arr[2][5]
        cov_low = arr[3][4]
        cov_mean = arr[3][1]
        cov_high = arr[3][5]
        row = [
                sequence_length, midpoint,
                mean_low, mean_mean, mean_high,
                var_low, var_mean, var_high,
                cov_low, cov_mean, cov_high]
        data_arr.append(row)
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table string
    table_string = RUtil.get_table_string(data_arr, g_headers)
    # get the scripts
    scripts = beasttut.get_ggplot2_scripts(
            nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, scripts

def get_loganalysis_array(start_pos, stop_pos, nsamples, header_seq_pairs):
    # get the xml contents
    xml_string = beasttut.get_xml_string(
            start_pos, stop_pos, nsamples, 'myjob.log', header_seq_pairs)
    # prepare the base path for the beast analysis
    basepath = beast.prepare()
    with open(os.path.join(basepath, 'myjob.xml'), 'w') as fout:
        fout.write(xml_string + '\n')
    beast.run_beast(basepath, g_beast_jar_path)
    beast.run_loganalyser(basepath)
    # read the analysis
    with open(os.path.join(basepath, 'myjob-loganalyser.txt')) as fin:
        analysis_text = fin.read()
    # parse the analysis
    arr = beast.loganalyser_to_array(analysis_text)
    return arr

def get_response_content(fs):
    start_pos = 1
    stop_pos = fs.ncolumns
    nsamples = fs.nsamples
    header_seq_pairs = beasttut.get_456_col_permuted_header_seq_pairs()
    arr = get_loganalysis_array(
            start_pos, stop_pos, nsamples, header_seq_pairs)
    s = '\n'.join('\t'.join(str(x) for x in row) for row in arr) + '\n'
    return s

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
    f.info('(local) permute columns of the alignment')
    header_seq_pairs = beasttut.get_456_col_permuted_header_seq_pairs()
    f.info('(local) run BEAST serially locally and build the R stuff')
    table_string, scripts = get_table_string_and_scripts(
            g_start_stop_pairs, args.nsamples, header_seq_pairs)
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
    """
    parser.add_argument('--run', choices=('serial', 'parallel', 'hpc'),
            default='serial',
            help='execution model')
    """
    parser.add_argument('-v', '--verbose',
            action='store_true',
            help='show more info')
    main(parser.parse_args())

