"""
Use BEAST to analyze multiple partitions of an alignment.

The idea is that different intervals will have different
estimates for the coefficient of variation
of substitution rates among branches.
The xml used by the web interface is a simplification
of the example BEAST primates.xml file.
The command line interface allows any xml file.
"""

from StringIO import StringIO
import argparse
import os
import sys
import logging

from lxml import etree
import gmpy

import Form
import FormOut
import beast
import beasttut
import beasttiling
import RUtil
import moretypes


g_beast_jar_path = os.path.expanduser('~/BEASTv1.7.1/lib/beast.jar')

g_ncols_max = 456

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
    # The web interface uses a fixed xml
    # and looks at only a single column interval.
    return [
            Form.IntegerInterval(
                'start_pos', 'stop_pos', 'alignment interval',
                1, g_ncols_max, low=1, high=g_ncols_max),
            Form.Integer('nsamples', 'mcmc chain steps',
                8000, low=80, high=8000)
            ]

def get_form_out():
    return FormOut.Report('summary')

def get_response_content(fs):
    header_sequence_pairs = beasttut.get_header_seq_pairs()
    interval_xml_data = beasttut.get_xml_string(
            fs.start_pos, fs.stop_pos, fs.nsamples, 'myjob.log',
            header_sequence_pairs)
    arr = get_loganalysis_array(interval_xml_data)
    s = '\n'.join('\t'.join(str(x) for x in row) for row in arr) + '\n'
    return s

def get_table_strings_and_scripts(
        xmldata, alignment_id, start_stop_pairs,
        nsamples):
    """
    Command-line only.
    @param xmldata: xml data already adjusted for nsamples and log filename
    @param alignment_id: xml element id
    @param start_stop_pairs: alignment interval bounds
    @param nsamples: an extra parameter for script generation
    @return: short table string, long table string, scripts (for short table)
    """
    # init the array for the full R table
    full_data_arr = []
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_pos, stop_pos in start_stop_pairs:
        sequence_length = stop_pos - start_pos + 1
        midpoint = (start_pos + stop_pos) / 2.0
        interval_xml_data = beast.set_alignment_interval(
                xmldata, alignment_id, start_pos, stop_pos)
        row_labels, col_labels, arr = get_loganalysis_labeled_array(
                interval_xml_data)
        stat_name_to_row_index = dict(
                (x, i) for i, x in enumerate(row_labels))
        summary_name_to_col_index = dict(
                (x, i) for i, x in enumerate(col_labels))
        # define row indices of interest
        mean_row_index = stat_name_to_row_index['meanRate']
        var_row_index = stat_name_to_row_index['coefficientOfVariation']
        cov_row_index = stat_name_to_row_index['covariance']
        # define column indices of interest
        mean_col_index = summary_name_to_col_index['mean']
        low_col_index = summary_name_to_col_index['hpdLower']
        high_col_index = summary_name_to_col_index['hpdUpper']
        row = [
                sequence_length,
                midpoint,
                arr[mean_row_index][low_col_index],
                arr[mean_row_index][mean_col_index],
                arr[mean_row_index][high_col_index],
                arr[var_row_index][low_col_index],
                arr[var_row_index][mean_col_index],
                arr[var_row_index][high_col_index],
                arr[cov_row_index][low_col_index],
                arr[cov_row_index][mean_col_index],
                arr[cov_row_index][high_col_index],
                ]
        data_arr.append(row)
        # add rows to the full data array
        for row_index, row_label in enumerate(row_labels):
            for col_index, col_label in enumerate(col_labels):
                row = [
                        sequence_length,
                        midpoint,
                        '"' + row_label + '"',
                        '"' + col_label + '"',
                        arr[row_index][col_index],
                        ]
                full_data_arr.append(row)
        # add entries to some utility arrays
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table strings
    table_string = RUtil.get_table_string(data_arr, g_headers)
    full_table_string = RUtil.get_table_string(
            full_data_arr,
            [
                'sequence.length',
                'midpoint',
                'statistic.name',
                'posterior.analysis',
                'value'
                ],
            force_float=False,
            )
    # get the scripts
    scripts = beasttut.get_ggplot2_scripts(
            nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, full_table_string, scripts

def get_loganalysis_text(xmldata):
    # prepare the base path for the beast analysis
    basepath = beast.prepare()
    with open(os.path.join(basepath, 'myjob.xml'), 'w') as fout:
        fout.write(xmldata + '\n')
    beast.run_beast(basepath, g_beast_jar_path)
    beast.run_loganalyser(basepath)
    # read the analysis
    with open(os.path.join(basepath, 'myjob-loganalyser.txt')) as fin:
        analysis_text = fin.read()
    return analysis_text

def get_loganalysis_array(xmldata):
    analysis_text = get_loganalysis_text(xmldata)
    return beast.loganalyser_to_array(analysis_text)

def get_loganalysis_labeled_array(xmldata):
    analysis_text = get_loganalysis_text(xmldata)
    return beast.loganalyser_to_labeled_array(analysis_text)

def main(args):
    # check args
    if gmpy.popcount(args.ntiles) != 1:
        raise ValueError('the number of tiles should be a power of two')
    # set up the logger
    f = logging.getLogger('toplevel.logger')
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter('%(message)s %(asctime)s'))
    f.addHandler(h)
    if args.verbose:
        f.setLevel(logging.DEBUG)
    else:
        f.setLevel(logging.WARNING)
    f.info('(local) read the xml contents')
    if args.infile is None:
        xmldata = sys.stdin.read()
    else:
        with open(args.infile) as fin:
            xmldata = fin.read()
    f.info('(local) modify the log filename and chain length xml contents')
    xmldata = beast.set_nsamples(xmldata, args.mcmc_id, args.nsamples)
    xmldata = beast.set_log_filename(xmldata, args.log_id, args.log_filename)
    xmldata = beast.set_log_logevery(xmldata, args.log_id, args.log_logevery)
    f.info('(local) define the hierarchically nested intervals')
    start_stop_pairs = tuple(
            (a+1,b) for a, b in beasttiling.gen_hierarchical_slices(
                args.tile_width, args.offset, args.tile_width * args.ntiles))
    f.info('(local) run BEAST serially locally and build the R stuff')
    table_string, full_table_string, scripts = get_table_strings_and_scripts(
            xmldata, args.alignment_id, start_stop_pairs, args.nsamples)
    if args.full_table_out:
        f.info('(local) create the verbose R table')
        with open(args.full_table_out, 'w') as fout:
            fout.write(full_table_string)
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
    # The command line mode looks at a hierarchical tiling within an interval.
    # Request a tile width and a number of tiles and an offset.
    # The mcmc and alignment element ids will need to be requested
    # before they can be identified within the xml;
    # this is so that the right ones can be modified.
    parser = argparse.ArgumentParser(description=__doc__)
    #
    # define input and output and verbosity
    parser.add_argument('-i', '--infile',
            help='input xml file')
    parser.add_argument('-o', '--outfile',
            default='beast-analysis.pdf',
            help='write this pdf file')
    parser.add_argument('-v', '--verbose',
            action='store_true',
            help='show more info')
    #
    # optionally write a more verbose R table
    parser.add_argument('--full_table_out',
            help='write a verbose R table to this location')
    #
    # define the mcmc id and the number of samples
    parser.add_argument('--mcmc_id',
            default='mcmc',
            help='xml mcmc element id')
    parser.add_argument('--nsamples',
            default=8000, type=moretypes.pos_int,
            help='let the BEAST MCMC generate this many samples')
    #
    # define the log id and the filename
    parser.add_argument('--log_id',
            default='fileLog',
            help='xml log element id')
    parser.add_argument('--log_filename',
            default='myjob.log',
            help='name of the log file to be read by loganalyser')
    parser.add_argument('--log_logevery',
            default=1, type=moretypes.pos_int,
            help='log the statistics once per this many mcmc samples')
    #
    # define the alignment id and the interval
    parser.add_argument('--alignment_id',
            default='alignment',
            help='xml alignment element id')
    parser.add_argument('--tile_width',
            default=50, type=moretypes.pos_int,
            help='smallest number of alignment columns per analysis')
    parser.add_argument('--ntiles',
            default=8, type=moretypes.pos_int,
            help='analyze an interval this many tiles wide')
    parser.add_argument('--offset',
            default=0, type=moretypes.nonneg_int,
            help='use this offset from the first alignment column')
    #
    # run the analysis
    main(parser.parse_args())

