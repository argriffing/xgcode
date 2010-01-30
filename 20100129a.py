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


import StringIO
import argparse
import os

from SnippetUtil import HandlingError
import Form


g_sample_lines = [
        'YHet 3261 T C/C 2 A 0 C 1 G 0 T 1 15 15 50',
        'YHet 4197 G C/C 2 A 0 C 1 G 1 T 0 2 2 60',
        'YHet 4573 C T/T 3 A 0 C 0 G 0 T 3 36 36 51',
        'YHet 5490 G A/A 7 A 4 C 0 G 3 T 0 2 2 0',
        '2L 5091 T C/T 16 A 0 C 7 G 0 T 9 27 27 28',
        '2L 5092 C C/T 16 A 0 C 9 G 0 T 7 30 78 28',
        '2L 5095 T A/T 17 A 8 C 0 G 0 T 9 31 82 31']


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
    out = StringIO.StringIO()
    print >> out, 'hello world'
    return [('Content-Type', 'text/plain')], out.getvalue().strip()

def gen_stripped_lines(line_source):
    for line in line_source:
        line = line.strip()
        if line:
            yield line

def line_to_tuple(line):
    """
    Parse a line and do error checking.
    @param line: a stripped line of input
    @return: a tuple with proper types
    """
    values = line.split()
    if len(values) != 16:
        raise Exception('expected 16 values per line')
    e = Exception('literal A, C, G, T letters were not found where expected')
    if values[5] != 'A' or values[7] != 'C':
        raise e
    if values[9] != 'G' or values[11] != 'T':
        raise e
    typed_values = [
            values[0], int(values[1]), values[2],
            values[3], int(values[4]),
            values[5], int(values[6]), values[7], int(values[8]),
            values[9], int(values[10]), values[11], int(values[12]),
            int(values[13]), int(values[14]), int(values[15])]
    return typed_values

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
    # scan the input file for chromosome names and formatting
    name_to_last_pos = {}
    with open(input_filename) as fin:
        for line in gen_stripped_lines(fin):
            values = line_to_tuple(line)
            name, pos = values[0], values[1]
            last_pos = name_to_last_pos.get(name, None)
            msg = 'expected strictly increasing positions per chromosome'
            if last_pos is not None:
                if last_pos >= pos:
                    raise Exception(msg)
            name_to_last_pos[name] = pos
    chromosome_names = list(sorted(name_to_last_pos.keys()))
    # define the full output file paths
    ch_paths = []
    for name in chromosome_names:
        a = output_directory
        b = os.path.basename(input_filename)
        c = '_' + name
        ch_path = os.path.join(a, b+c)
        ch_paths.append(ch_path)
    # check for existence of output
    if not force:
        for p in ch_paths:
            if os.path.exists(p):
                raise Exception('output already exists: ' + p)
    # create the files open for writing
    ch_files = []
    for p in ch_paths:
        ch_files.append(open(p, 'wt'))
    # write the headers
    header = '\t'.join(['position', 'A', 'C', 'G', 'T'])
    for f in ch_files:
        f.write(header + '\n')
    # write the lines
    name_to_file = dict(zip(chromosome_names, ch_files))
    with open(input_filename) as fin:
        for line_in in gen_stripped_lines(fin):
            V = line_to_tuple(line_in)
            name, pos = V[0], V[1]
            A, C, G, T = V[6], V[8], V[10], V[12]
            f = name_to_file[name]
            line_out = '\t'.join(str(x) for x in (pos, A, C, G, T))
            f.write(line_out + '\n')
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
    args = parser.parse_args()
    main(args)
