"""Convert a filtered pileup file to a set of files.

Each output file has information from a single chromosome.
Input columns are
{chromosome name, chromosome position, reference base,
calls for the two alleles, coverage, literal A, A count, literal C,
C count, literal G, G count, literal T, T count, first quality score,
second quality score, third quality score}.
Output includes a header line followed by data lines.
The output columns are
{position, A count, C count, G count, T count}.
The output filenames include the chromosome name.
"""


from StringIO import StringIO
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Progress
import DGRP
import iterutils
import const

g_sample_data = const.read('20100730c')

g_header = '\t'.join(['position', 'A', 'C', 'G', 'T'])


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('data_in', 'filtered pileup file', g_sample_data)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # quickly skim the lines to get some info
    fin = StringIO(fs.data_in)
    skimmer = DGRP.ChromoSkimmer()
    for chromo_name in skimmer.skim(gen_untyped_rows(fin)):
        pass
    chromo_names = skimmer.name_list
    nlines = skimmer.linecount
    # check formatting and monotonicity
    fin = StringIO(fs.data_in)
    for i in DGRP.check_chromo_monotonicity(gen_typed_rows(fin)):
        pass
    # begin writing
    out = StringIO()
    print >> out, 'writing the first of', len(chromo_names), 'chromosomes:'
    print >> out
    # write only the first chromosome
    fin = StringIO(fs.data_in)
    print >> out, g_header
    for row in gen_typed_rows(fin):
        name = row[0]
        if name == chromo_names[0]:
            print >> out, '\t'.join(str(x) for x in convert_row(row))
    return [('Content-Type', 'text/plain')], out.getvalue().strip()

def line_to_row(line):
    """
    Parse a line and do error checking.
    @param line: a stripped line of input
    @return: a tuple with proper types
    """
    values = line.split()
    if len(values) != 16:
        raise Exception('expected 16 values per line')
    msg = 'literal A, C, G, T letters were not found where expected'
    if values[5] != 'A' or values[7] != 'C':
        raise Exception(msg)
    if values[9] != 'G' or values[11] != 'T':
        raise Exception(msg)
    typed_values = [
            values[0], int(values[1]), values[2],
            values[3], int(values[4]),
            values[5], int(values[6]), values[7], int(values[8]),
            values[9], int(values[10]), values[11], int(values[12]),
            int(values[13]), int(values[14]), int(values[15])]
    return typed_values

def gen_typed_rows(fin):
    for line in iterutils.stripped_lines(fin):
        yield line_to_row(line)

def gen_untyped_rows(fin):
    for line in iterutils.stripped_lines(fin):
        yield line.split()

def convert_row(row):
    """
    @param row: a sequence of values
    @return: a subset of the values
    """
    name, pos = row[0], row[1]
    A, C, G, T = row[6], row[8], row[10], row[12]
    return (pos, A, C, G, T)

def main(args):
    """
    @param args: positional and flaglike arguments
    """
    # read the arguments
    input_filename = os.path.abspath(os.path.expanduser(args.infile))
    output_directory = os.path.abspath(os.path.expanduser(args.outdir))
    force = args.force
    # make sure that the output directory exists
    if not os.path.isdir(output_directory):
        if force:
            os.makedirs(output_directory)
    if not os.path.isdir(output_directory):
        msg = 'output directory does not exist: ' + output_directory
        raise Exception(msg)
    # scan the input file for chromosome names
    ch_paths = []
    skimmer = DGRP.ChromoSkimmer()
    with open(input_filename) as fin:
        for chromo_name in skimmer.skim(gen_untyped_rows(fin)):
            output_filename = args.out_prefix + chromo_name + args.out_suffix
            ch_path = os.path.join(output_directory, output_filename)
            ch_paths.append(ch_path)
            if not force:
                if os.path.exists(ch_path):
                    raise Exception('output already exists: ' + ch_path)
    chromo_names = skimmer.name_list
    nlines = skimmer.linecount
    # start the progress bar
    nticks = 2*nlines
    pbar = Progress.Bar(nticks)
    # scan the input file for correct types and for monotonicity
    with open(input_filename) as fin:
        for i in DGRP.check_chromo_monotonicity(gen_typed_rows(fin)):
            pbar.increment()
    # create the files open for writing
    ch_files = []
    for p in ch_paths:
        ch_files.append(open(p, 'wt'))
    # write the headers
    if not args.noheader:
        for f in ch_files:
            f.write(g_header + '\n')
    # write the lines
    name_to_file = dict(zip(chromo_names, ch_files))
    with open(input_filename) as fin:
        for row in gen_typed_rows(fin):
            name = row[0]
            row_out = convert_row(row)
            f = name_to_file[name]
            line_out = '\t'.join(str(x) for x in row_out)
            f.write(line_out + '\n')
            pbar.increment()
    # close the files
    for f in ch_files:
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--outdir', default=os.getcwd(),
            help='write the chromosome files to this directory')
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--noheader', action='store_true',
            help='omit the header row from the output'),
    parser.add_argument('--out_prefix', default='chromosome.',
            help='prefix added to the chromosome name in the output filename')
    parser.add_argument('--out_suffix', default='.txt',
            help='suffix added to the chromosome name in the output filename')
    args = parser.parse_args()
    main(args)
