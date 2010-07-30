"""Convert a filtered pileup file to files with a minimal observation per line.

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


from StringIO import StringIO
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import Progress
import DGRP
import ambignt
import iterutils
import iterfiller
import const

g_sample_data = const.read('20100730b')

class Scanner:
    """
    Go through a filtered pileup file and check and save chromosome info.
    """
    def __init__(self, low, high, fill, errlow, errhigh):
        """
        @param first: an integer or None or 'drosophila'
        @param last: an integer or None or 'drosophila'
        @param fill: True to fill gaps with default values
        @param errlow: True to flag low positions
        @param errhigh: True to flag high positions
        """
        self.low = low
        self.high = high
        self.errlow = errlow
        self.errhigh = errhigh
        self.fill = fill
        self.name_to_counter = {}
        self.name_to_generator = {}

    def get_npositions(self):
        """
        @return: the total number of lines to write
        """
        return sum(c.npositions for c in self.name_to_counter.values())

    def scan(self, fin):
        """
        Save chromosome info and check for errors.
        Yield chromosome names as they are encountered.
        @param fin: a file open for reading
        """
        name_to_drosophila_length = dict(DGRP.g_chromosome_length_pairs)
        last_name = None
        for row in gen_typed_rows(fin):
            name, position = row[0], row[1]
            # assert that chromosomes are contiguous
            if name != last_name:
                if name in self.name_to_counter:
                    msg = 'chromosome ' + name + ' should be contiguous'
                    raise Exception(msg)
            # create info for a new chromosome if necessary
            if name not in self.name_to_counter:
                # define the high value
                if self.high == 'drosophila':
                    high = name_to_drosophila_length.get(name, None)
                    if high is None:
                        raise Exception('invalid fly chromosome: ' + name)
                else:
                    high = self.high
                # define the low value
                if self.low == 'drosophila':
                    low = 1
                else:
                    low = self.low
                # create the counter and the generator
                fc = iterfiller.FillerCounter(low, high, self.fill,
                        self.errlow, self.errhigh)
                fg = iterfiller.FillerGenerator(low, high, self.fill,
                        self.errlow, self.errhigh, None)
                self.name_to_counter[name] = fc
                self.name_to_generator[name] = fg
                # update the current chromosome name
                last_name = name
                # yield the new name so the file can be checked
                yield name
            # use the counter to add the position
            self.name_to_counter[name].fill(position)
        for c in self.name_to_counter.values():
            c.finish()

    def gen_named_observations(self, fin):
        """
        Yield (chrom_name, observation) pairs
        @param fin: a file open for reading
        """
        default_value = None
        # Process each row of the input file,
        # yielding after each written line.
        for row in gen_typed_rows(fin):
            name, position = row[0], row[1]
            value = DGRP.filtered_pileup_typed_to_obs(row)
            fg = self.name_to_generator[name]
            for obs in fg.fill(position, value):
                yield name, obs
        for name, fg in self.name_to_generator.items():
            for obs in fg.finish():
                yield name, obs

    def gen_named_lines(self, fin):
        """
        Yield (chrom_name, observation_line) pairs
        @param fin: a file open for reading
        """
        default_obs = (0, 0, 0, 0)
        default_line = '\t'.join(str(x) for x in default_obs)
        # yield chromosome names and observation lines
        for name, obs in self.gen_named_observations(fin):
            if obs is None:
                line = default_line
            else:
                line = '\t'.join(str(x) for x in obs)
            yield (name, line)


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('data_in', 'filtered pileup file', g_sample_data),
            Form.RadioGroup('low_info', 'low position', [
                Form.RadioItem('low_0', '0'),
                Form.RadioItem('low_1', '1', True),
                Form.RadioItem('low_none', 'none')]),
            Form.RadioGroup('high_info', 'high position', [
                Form.RadioItem('high_1000', '1000'),
                Form.RadioItem('high_none', 'none', True)]),
            Form.CheckGroup('options', 'output options', [
                Form.CheckItem('errlow', 'err on out of bounds low', True),
                Form.CheckItem('errhigh', 'err on out of bounds high', True),
                Form.CheckItem('fill', 'fill with default values', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # unpack the first and last requested positions
    low = {'low_0':0, 'low_1':1, 'low_none':None}[fs.low_info]
    high = {'high_1000':1000, 'high_none':None}[fs.high_info]
    # create the scanner object which will be used for two passes
    scanner = Scanner(low, high, fs.fill, fs.errlow, fs.errhigh)
    # Do the first pass; check for errors and gather chromosome info.
    names = set()
    fin = StringIO(fs.data_in)
    for name in scanner.scan(fin):
        names.add(name)
    names = list(sorted(names))
    # See if the number of lines to be written is appropriate.
    npos = scanner.get_npositions()
    if npos > 2000:
        msg_a = 'attempting to write too many lines: '
        msg_b = '%d lines in %d files.' % (npos, len(names))
        raise HandlingError(msg_a + msg_b)
    # Do the second pass; write the response for only the first chromosome
    out = StringIO()
    print >> out, 'writing the first of', len(names), 'chromosomes:'
    print >> out
    fin = StringIO(fs.data_in)
    for name, line in scanner.gen_named_lines(fin):
        if name == names[0]:
            print >> out, line
    return [('Content-Type', 'text/plain')], out.getvalue().strip()

def gen_typed_rows(fin):
    for line in fin:
        srow = line.split()
        if srow:
            yield DGRP.filtered_pileup_row_to_typed(srow)

def main(args):
    # read the arguments
    input_filename = os.path.abspath(os.path.expanduser(args.infile))
    output_directory = os.path.abspath(os.path.expanduser(args.outdir))
    force = args.force
    low, high = args.low, args.high
    errlow, errhigh = args.errlow, args.errhigh
    # make sure that the output directory exists
    if not os.path.isdir(output_directory):
        if force:
            os.makedirs(output_directory)
    if not os.path.isdir(output_directory):
        msg = 'output directory does not exist: ' + output_directory
        raise Exception(msg)
    # create the scanner object which will be used for two passes
    scanner = Scanner(low, high, args.fill, errlow, errhigh)
    # Do the first pass,
    # checking for errors and gathering info about the chromosomes.
    name_to_path = {}
    with open(input_filename) as fin:
        for name in scanner.scan(fin):
            output_filename = args.out_prefix + name + args.out_suffix
            fpath = os.path.join(output_directory, output_filename)
            name_to_path[name] = fpath
            if not args.force:
                if os.path.exists(fpath):
                    raise Exception('output file already exists: ' + fpath)
    nticks = scanner.get_npositions()
    pbar = Progress.Bar(nticks)
    # open the files for writing
    name_to_fout = {}
    for name, fpath in name_to_path.items():
        name_to_fout[name] = open(fpath, 'wt')
    # Do the second pass,
    # writing the files and updating the progress bar.
    with open(input_filename) as fin:
        for name, line in scanner.gen_named_lines(fin):
            name_to_fout[name].write(line + '\n')
            pbar.increment()
    # close the files
    for fout in name_to_fout.values():
        fout.close()

def drosophila_position(value):
    """
    This is a argparse compatible type.
    """
    try:
        v = int(value)
    except ValueError, e:
        v = None
    if v is None:
        if value == 'none':
            return None
        elif value == 'drosophila':
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
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--fill', action='store_true',
            help='fill missing positions with default values')
    parser.add_argument('--errlow', action='store_true',
            help='assert that no position is below the lower bound')
    parser.add_argument('--errhigh', action='store_true',
            help='assert that no position is above the upper bound')
    parser.add_argument('--low', default='drosophila',
            type=drosophila_position,
            metavar='{<int>, drosophila, none}',
            help='the first position in a chromosome')
    parser.add_argument('--high', default='drosophila',
            type=drosophila_position,
            metavar='{<int>, drosophila, none}',
            help='the last position in a chromosome')
    parser.add_argument('--outdir', default=os.getcwd(),
            help='write the chromosome files to this directory')
    parser.add_argument('--out_prefix', default='chromosome.',
            help='prefix added to the chromosome name in the output filename')
    parser.add_argument('--out_suffix', default='.txt',
            help='suffix added to the chromosome name in the output filename')
    main(parser.parse_args())
