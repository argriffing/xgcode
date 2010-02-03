"""Convert a filtered pileup file to a set of files.

Each output file has information from a single chromosome.
Input columns are
{chromosome name, chromosome position, reference base,
calls for the two alleles, coverage, literal A, A count, literal C,
C count, literal G, G count, literal T, T count, first quality score,
second quality score, third quality score}.
The chromosome position numbers may be indexed from 0 or from 1.
Each line of output corresponds to a chromosome position,
starting from the first position and proceding consecutively,
padded with zeros.
The total chromosome length may be provided.
Output lines are four whitespace separated integers.
The first integer is the number of reads at the position
aligned to the reference nucleotide.
The remaining three integers are the numbers of reads
aligned to the non-reference nucleotides.
These three nucleotides are sorted in decreasing order.
No header line is provided.
"""


import StringIO
import argparse
import os

from SnippetUtil import HandlingError
import Form
import Progress
import DGRP
import Util


g_sample_lines = [
        'chrI 61 T C/C 2 A 0 C 1 G 0 T 1 15 15 50',
        'chrI 67 G C/C 2 A 0 C 1 G 1 T 0 2 2 60',
        'chrI 73 C T/T 3 A 0 C 0 G 0 T 3 36 36 51',
        'chrI 90 G A/A 7 A 4 C 0 G 3 T 0 2 2 0',
        'chrQ 91 T C/T 16 A 0 C 7 G 0 T 9 27 27 28',
        'chrQ 92 C C/T 16 A 0 C 9 G 0 T 7 30 78 28',
        'chrQ 95 T A/T 17 A 8 C 0 G 0 T 9 31 82 31']



def get_form():
    """
    @return: the body of a form
    """
    sample_data = '\n'.join(g_sample_lines)
    form_objects = [
            Form.MultiLine('data_in', 'filtered pileup file', sample_data)]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # quickly skim the lines to get some info
    fin = StringIO.StringIO(fs.data_in)
    skimmer = DGRP.ChromoSkimmer()
    for chromo_name in skimmer.skim(gen_untyped_rows(fin)):
        pass
    chromo_names = skimmer.name_list
    nlines = skimmer.linecount
    # check formatting and monotonicity
    fin = StringIO.StringIO(fs.data_in)
    for i in DGRP.check_chromo_monotonicity(gen_typed_rows(fin)):
        pass
    # begin writing
    out = StringIO.StringIO()
    print >> out, 'writing the first of', len(chromo_names), 'chromosomes:'
    print >> out
    # write only the first chromosome
    fin = StringIO.StringIO(fs.data_in)
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
    for line in Util.stripped_lines(fin):
        yield line_to_row(line)

def gen_untyped_rows(fin):
    for line in Util.stripped_lines(fin):
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
    nticks = args.length
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

def first_position(value):
    try:
        v = int(value)
    except ValueError, e:
        v = None
    if v is None:
        if value in ('min', 'drosophila'):
            return value
        else:
            raise TypeError()
    else:
        if v < 0:
            raise TypeError()
        else:
            return v

def last_position(value):
    try:
        v = int(value)
    except ValueError, e:
        v = None
    if v is None:
        if value in ('max', 'drosophila'):
            return value
        else:
            raise TypeError()
    else:
        if v < 0:
            raise TypeError()
        else:
            return v

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--outdir', default=os.getcwd(),
            help='write the chromosome files to this directory')
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--out_prefix', default='chromosome.',
            help='prefix added to the chromosome name in the output filename')
    parser.add_argument('--out_suffix', default='.txt',
            help='suffix added to the chromosome name in the output filename')
    parser.add_argument('--first', default='drosophila', type=first_position,
            metavar='{<int>, min, drosophila}',
            help='the first position in a chromosome'),
    parser.add_argument('--last', default='drosophila', type=last_position,
            metavar='{<int>, max, drosophila}',
            help='the last position in a chromosome'),
    args = parser.parse_args()
    main(args)
