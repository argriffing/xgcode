"""
A module for interfacing with BEAST.
"""

# FIXME UNUSED

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


g_remote_beast_sh_path = os.path.join(
        '/brc_share/brc/argriffi/packages/BEASTv1.7.1',
        'bin/beast')
g_remote_beast_jar_path = os.path.join(
        '/brc_share/brc/argriffi/packages/BEASTv1.7.1',
        'lib/beast.jar')
g_remote_java_path = '/usr/local/apps/java/jre1.6.0_31/bin/java'



class RemoteBeast(hpcutil.RemoteBrc):
    def __init__(self, start_stop_pairs, nsamples):
        hpcutil.RemoteBrc.__init__(self)
        self.start_stop_pairs = start_stop_pairs
        self.nsamples = nsamples
        self.local_log_paths = None
    def preprocess(self):
        # for each (start_pos, stop_pos) pair
        # define a log filename and
        # create a beast xml and a bsub script
        self.local_log_paths = []
        for i, start_stop_pair in enumerate(self.start_stop_pairs):
            start_pos, stop_pos = start_stop_pair
            # define the log file paths
            log_name = 'beast-%d-%d.log' % start_stop_pair
            local_log_path = os.path.join(
                    self.local.get_out(), log_name)
            remote_log_path = os.path.join(
                    self.remote.get_out(), log_name)
            self.local_log_paths.append(local_log_path)
            # define the xml file paths
            xml_name = 'beast-%d-%d.xml' % start_stop_pair
            local_xml_path = os.path.join(
                    self.local.get_in_contents(), xml_name)
            remote_xml_path = os.path.join(
                    self.remote.get_in_contents(), xml_name)
            # define the local bsub path
            bsub_name = 'beast-%d-%d.bsub' % start_stop_pair
            local_bsub_path = os.path.join(
                    self.local.get_in_bsubs(), bsub_name)
            # create the xml file
            with open(local_xml_path, 'w') as fout:
                print >> fout, beasttut.get_xml_string(
                        start_pos, stop_pos, self.nsamples, remote_log_path)
            # create the bsub file
            with open(local_bsub_path, 'w') as fout:
                stdout_path = os.path.join(self.remote.get_out(),
                        'out.beast-%d-%d' % start_stop_pair)
                stderr_path = os.path.join(self.remote.get_out(),
                        'err.beast-%d-%d' % start_stop_pair)
                # name the job
                print >> fout, '#BSUB -J %d-%d' % start_stop_pair
                # suggest the brc queue
                print >> fout, '#BSUB -q brc'
                # redirect stdout
                print >> fout, '#BSUB -o', stdout_path
                # redirect stderr
                print >> fout, '#BSUB -e', stderr_path
                # run BEAST using a java path suggested by Gary Howell
                print >> fout, g_remote_java_path, '-jar'
                print >> fout, g_remote_beast_jar_path, remote_xml_path

def get_table_string_and_scripts(start_stop_pairs, nsamples):
    """
    Command-line only.
    """
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_pos, stop_pos in start_stop_pairs:
        sequence_length = stop_pos - start_pos + 1
        means, variations, covs = get_value_lists(
                start_pos, stop_pos, nsamples)
        midpoint = (start_pos + stop_pos) / 2.0
        row = [sequence_length, midpoint]
        for values in means, variations, covs:
            corr_info = mcmc.Correlation()
            corr_info.analyze(values)
            hpd_low, hpd_high = mcmc.get_hpd_interval(0.95, values)
            row.extend([hpd_low, corr_info.mean, hpd_high])
        data_arr.append(row)
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table string
    table_string = RUtil.get_table_string(data_arr, g_headers)
    # get the scripts
    scripts = get_ggplot2_scripts(nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, scripts

def forked_function(start_stop_n):
    """
    This function should accept and return as little data as possible.
    In particular do not return a huge multimegabyte nested list.
    @param start_stop_n: start_pos, stop_pos, nsamples
    @return: (corr_info, hpd_interval) for mean, variation, covariance
    """
    start_pos, stop_pos, nsamples = start_stop_n
    means, variations, covs = get_value_lists(start_pos, stop_pos, nsamples)
    post_pairs = []
    for values in means, variations, covs:
        corr_info = mcmc.Correlation()
        corr_info.analyze(values)
        hpd_interval = mcmc.get_hpd_interval(0.95, values)
        post_pairs.append((corr_info, hpd_interval))
    return post_pairs

def get_table_string_and_scripts_par(start_stop_pairs, nsamples):
    """
    Local command-line multi-process only.
    """
    # define the pool of processes corresponding to the number of cores
    mypool = Pool(processes=4)
    # do the multiprocessing
    start_stop_n_triples = [(a, b, nsamples) for a, b in start_stop_pairs]
    post_pairs_list = mypool.map(forked_function, start_stop_n_triples)
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_stop_pair, post_pairs in zip(start_stop_pairs, post_pairs_list):
        start_pos, stop_pos = start_stop_pair
        sequence_length = stop_pos - start_pos + 1
        midpoint = (start_pos + stop_pos) / 2.0
        row = [sequence_length, midpoint]
        for corr_info, hpd_interval in post_pairs:
            hpd_low, hpd_high = hpd_interval
            row.extend([hpd_low, corr_info.mean, hpd_high])
        data_arr.append(row)
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table string
    table_string = RUtil.get_table_string(data_arr, g_headers)
    # get the scripts
    scripts = get_ggplot2_scripts(nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, scripts

def get_table_string_and_scripts_from_logs(
        start_stop_pairs, log_paths, nsamples):
    """
    This is for analysis of remote execution.
    """
    # build the array for the R table
    data_arr = []
    sequence_lengths = []
    midpoints = []
    for start_stop_pair, log_path in zip(
            start_stop_pairs, log_paths):
        start_pos, stop_pos = start_stop_pair
        sequence_length = stop_pos - start_pos + 1
        means, variations, covs = read_log(log_path, nsamples)
        midpoint = (start_pos + stop_pos) / 2.0
        row = [sequence_length, midpoint]
        for values in means, variations, covs:
            corr_info = mcmc.Correlation()
            corr_info.analyze(values)
            hpd_low, hpd_high = mcmc.get_hpd_interval(0.95, values)
            row.extend([hpd_low, corr_info.mean, hpd_high])
        data_arr.append(row)
        sequence_lengths.append(sequence_length)
        midpoints.append(midpoint)
    # build the table string
    table_string = RUtil.get_table_string(data_arr, g_headers)
    # get the scripts
    scripts = get_ggplot2_scripts(nsamples, sequence_lengths, midpoints)
    # return the table string and scripts
    return table_string, scripts


