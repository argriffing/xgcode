"""
A module for interfacing with the BEAST MCMC java program.

This assumes that BEAST is on the path and that BEAGLE has been installed.
1) prepare() the base path
2) write the "<basepath>/myjob.xml" file referencing "myjob.log"
3) run_beast(basepath)
4) run_loganalyser(basepath)
The module also has a couple of functions to modify the BEAST xml files.
Perhaps these should be moved into their own module.
"""

from StringIO import StringIO
import subprocess
import time
import os
import tempfile

from lxml import etree

import Util

g_loganalyser_headers = [
        'statistic', 'mean', 'stdErr', 'median',
        'hpdLower', 'hpdUpper', 'ESS',
        '50hpdLower', '50hpdUpper']

def to_float(s):
    """
    The loganalyser output has annoying infinity symbols.
    I think it is ascii 236.
    """
    loganalyser_infinity = chr(226) + chr(136) + chr(158)
    python_infinity = 'inf'
    s = s.replace(loganalyser_infinity, python_infinity)
    try:
        value = float(s)
    except ValueError as e:
        raise ValueError(str(e) + ' : ' + ' '.join(str(ord(x)) for x in s))
    return value

def set_log_logevery(xmldata, log_id, logevery):
    """
    Sometimes the log id is fileLog.
    @param xmldata: xml file contents
    @param log_id: id of the log xml element
    @param logevery: log the statistics once every this many mcmc steps
    @return: new xml file contents
    """
    # read the xml tree
    tree = etree.parse(StringIO(xmldata))
    # modify the logging frequency
    for event, element in etree.iterwalk(tree, tag='log'):
        if element.get('id') == log_id:
            element.set('logEvery', str(logevery))
    # write the xml tree
    return etree.tostring(tree)

def set_log_filename(xmldata, log_id, filename):
    """
    Sometimes the log id is fileLog.
    @param xmldata: xml file contents
    @param log_id: id of the log xml element
    @param filename: name of the logfile
    @return: new xml file contents
    """
    # read the xml tree
    tree = etree.parse(StringIO(xmldata))
    # modify the log filename
    for event, element in etree.iterwalk(tree, tag='log'):
        if element.get('id') == log_id:
            element.set('fileName', str(filename))
    # write the xml tree
    return etree.tostring(tree)

def set_nsamples(xmldata, mcmc_id, nsamples):
    """
    @param xmldata: xml file contents
    @param mcmc_id: id of the mcmc xml element
    @param nsamples: target number of samples
    @return: new xml file contents
    """
    # read the xml tree
    tree = etree.parse(StringIO(xmldata))
    # modify the number of mcmc steps
    for event, element in etree.iterwalk(tree, tag='mcmc'):
        if element.get('id') == mcmc_id:
            element.set('chainLength', str(nsamples))
    # write the xml tree
    return etree.tostring(tree)

def _modify_taxon_sequence(taxon_element, start_pos, stop_pos):
    """
    This is a helper function for xml modification.
    """
    sequence = taxon_element.tail.strip()
    taxon_element.tail = sequence[start_pos-1 : stop_pos]

def set_alignment_interval(xmldata, alignment_id, start_pos, stop_pos):
    """
    @param xmldata: xml file contents
    @param alignment_id: id of the alignment xml element
    @param start_pos: one-based start position
    @param stop_pos: one-based inclusive stop position
    @return: new xml file contents
    """
    # read the xml tree
    tree = etree.parse(StringIO(xmldata))
    # modify the sequences within the alignment
    for event, element in etree.iterwalk(tree, tag='alignment'):
        if element.get('id') == alignment_id:
            for seq_element in element:
                if seq_element.tag != 'sequence':
                    continue
                for taxon_element in seq_element:
                    if taxon_element.tag != 'taxon':
                        continue
                    _modify_taxon_sequence(
                            taxon_element, start_pos, stop_pos)
    # write the xml tree
    return etree.tostring(tree)

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



"""
burnIn   <= 800,   maxState  = 8000
statistic   mean    stdErr  median  hpdLower    hpdUpper    ESS 50hpdLower  50hpdUpper
meanRate    0.01    1.0107E-3   9.8709E-3   8.4519E-3   0.0118  7.3723  8.6733E-3   0.0101  *
coefficientOfVariation  0.589   0.0989  0.5851  0.4099  0.7753  18.0456 0.5556  0.6628  *
covariance  -0.2115 0.113   -0.2303 -0.4392 -0.0295 5.1801  -0.3593 -0.1941 *

 * WARNING: The results of this MCMC analysis may be invalid as 
             one or more statistics had very low effective sample sizes (ESS)
"""

class LoganalyserParsingError(Exception): pass

def loganalyser_to_array(lstr):
    """
    @param lstr: contents of a loganalyser file
    @return: an array whose first row is headers
    """
    arr = []
    lines = [line.strip() for line in lstr.splitlines()]
    if not lines[0].startswith('burnIn'):
        raise LoganalyserParsingError(
                'expected the first line to start with burnIn')
    second_line_entries_observed = lines[1].split()
    if second_line_entries_observed != g_loganalyser_headers:
        raise LoganalyserParsingError(
                'expected the second line to have the headers; '
                'observed entries: %s '
                'expected entries: %s ' % (
                    second_line_entries_observed, g_loganalyser_headers))
    arr.append(g_loganalyser_headers)
    for line in lines[2:]:
        # stop if we reach a blank line
        if not line:
            break
        # add the row to the array
        row = line.split()[:len(g_loganalyser_headers)]
        try:
            row = [row[0]] + [to_float(x) for x in row[1:]]
        except ValueError as e:
            raise ValueError(str(e) + '\n' + lstr)
        arr.append(row)
    return arr

def loganalyser_to_labeled_array(lstr):
    """
    @param lstr: contents of a loganalyser file
    @return: row labels, column labels, array
    """
    arr_in = loganalyser_to_array(lstr)
    row_labels = zip(*arr_in)[0][1:]
    column_labels = arr_in[0][1:]
    arr_out = []
    for row_in in arr_in[1:]:
        row_out = row_in[1:]
        arr_out.append(row_out)
    return row_labels, column_labels, arr_out

