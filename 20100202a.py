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


import StringIO
import argparse
import os
import profile

from SnippetUtil import HandlingError
import Form
import Progress
import DGRP
import Util
import ambignt


g_sample_lines = [
        'chrI 61 T C/C 2 A 0 C 1 G 0 T 1 15 15 50',
        'chrI 67 G C/C 2 A 0 C 1 G 1 T 0 2 2 60',
        'chrI 73 C T/T 3 A 0 C 0 G 0 T 3 36 36 51',
        'chrI 90 G A/A 7 A 4 C 0 G 3 T 0 2 2 0',
        'chrQ 91 T C/T 16 A 0 C 7 G 0 T 9 27 27 28',
        'chrQ 92 C C/T 16 A 0 C 9 G 0 T 7 30 78 28',
        'chrQ 95 T A/T 17 A 8 C 0 G 0 T 9 31 82 31']


class Filler:
    """
    Generate values over a range of sequential positions.
    At some of the positions a value is available,
    but at other positions no value is available.
    This implementation is memory efficient because it uses iterators.
    """
    def __init__(self, low, high):
        """
        @param low: the position for which the first value is yielded
        @param high: the position for which the last value is yielded
        """
        self.low = low
        self.high = high
        self.prev = None

    def fill(self, position, value, default_value, finish=False):
        """
        Yield an informative value and maybe some uninformative ones.
        This function should be called repeatedly,
        and with strictly increasing positions.
        @param position: an available position
        @param value: the value at the position
        @param default_value: the value at missing positions
        @param finish: True if this is the last available value
        """
        # if the position is outside the range, just drop it
        if not (self.low <= position <= self.high):
            return
        # check monotonicity
        if self.prev is not None:
            if position <= self.prev:
                raise ValueError('positions should monotonically increase')
        # fill between the previous position and the current position
        for i in xrange(self.get_ngap(position)):
            yield default_value
        # yield the value at the current position
        yield value
        self.prev = position
        # possibly finish filling the range
        if finish:
            for i in xrange(self.get_nremaining()):
                yield default_value
            self.prev = self.high

    def get_ngap(self, position):
        if self.prev is None:
            return position - self.low
        else:
            return (position - self.prev) - 1

    def get_nremaining(self):
        return self.get_ngap(self.high) + 1



class ChromInfo:
    """
    Chromosome info.
    """
    def __init__(self, name):
        self.name = name
        self.low = None
        self.high = None
        self.nobserved = 0
    def add_position(self, position):
        self.nobserved += 1
        if self.low is None:
            self.low = position
            self.high = position
        else:
            if position <= self.high:
                raise ValueError('positions added out of order')
        self.low = min(self.low, position)
        self.high = max(self.high, position)

def get_requested_low(chrom, first):
    if first == 'drosophila':
        return 1
    elif first in ('min', 'ignore'):
        return chrom.low
    else:
        return first

def get_requested_high(chrom, last):
    if last in ('drosophila', 'truncated_drosophila'):
        return dict(DGRP.g_chromosome_length_pairs)[chrom.name]
    elif last in ('max', 'ignore'):
        return chrom.high
    else:
        return last

def get_requested_npositions(chrom, first, last):
    if first == 'ignore':
        return chrom.nobserved
    else:
        requested_low = get_requested_low(chrom, first)
        requested_high = get_requested_high(chrom, last)
        return (requested_high - requested_low) + 1


class Scanner:
    """
    Go through a filtered pileup file and check and save chromosome info.
    """
    def __init__(self, first, last):
        """
        @param first: an integer or 'min' or 'drosophila'
        @param last: an integer or 'max' or ['truncated_']'drosophila'
        """
        check_position_requests(first, last)
        self.first = first
        self.last = last
        self.name_to_chrom = {}

    def get_npositions(self):
        """
        @return: the total number of lines to write
        """
        chroms = self.name_to_chrom.values()
        first = self.first
        last = self.last
        return sum(get_requested_npositions(c, first, last) for c in chroms)

    def scan(self, fin):
        """
        Save chromosome info and check for errors.
        Yield chromosome names as they are encountered.
        @param fin: a file open for reading
        """
        last_name = None
        for row in gen_typed_rows(fin):
            name, position = row[0], row[1]
            # assert that chromosomes are contiguous
            if name != last_name:
                if name in self.name_to_chrom:
                    msg = 'chromosome ' + name + ' should be contiguous'
                    raise Exception(msg)
            # assert that drosophila-specific properties are correct
            if self.last in ('drosophila', 'truncated_drosophila'):
                name_to_length = dict(DGRP.g_chromosome_length_pairs)
                # assert that the chromosome is a valid Drosophila name
                if name not in name_to_length:
                    raise Exception('invalid Drosophila chromosome: ' + name)
                # assert that the position is not too high
                if self.last == 'drosophila':
                    if position > name_to_length[name]:
                        msg = 'position out of range: ' + str(position)
                        raise Exception(msg)
            # create info for a new chromosome if necessary
            if name not in self.name_to_chrom:
                self.name_to_chrom[name] = ChromInfo(name)
                last_name = name
                yield name
            # update the chromosome info with the position
            self.name_to_chrom[name].add_position(position)

    def gen_named_observations_unfilled(self, fin):
        """
        Yield (chrom_name, observation) pairs
        @param fin: a file open for reading
        """
        for row in gen_typed_rows(fin):
            name = row[0]
            chrom = self.name_to_chrom[name]
            observation = convert_row(row)
            yield name, observation

    def gen_named_observations_filled(self, fin):
        """
        Yield (chrom_name, observation) pairs
        @param fin: a file open for reading
        """
        default_value = None
        # Create a filler object for each chromosome.
        name_to_filler = {}
        for name, chrom in self.name_to_chrom.items():
            filler_low = get_requested_low(chrom, self.first)
            filler_high = get_requested_high(chrom, self.last)
            name_to_filler[name] = Filler(filler_low, filler_high)
        # Process each row of the input file,
        # yielding after each written line.
        for row in gen_typed_rows(fin):
            name, position = row[0], row[1]
            chrom = self.name_to_chrom[name]
            value = convert_row(row)
            filler = name_to_filler[name]
            finish = (position == chrom.high)
            for obs in filler.fill(position, value, default_value, finish):
                yield name, obs

    def gen_named_lines(self, fin):
        """
        Yield (chrom_name, observation_line) pairs
        @param fin: a file open for reading
        """
        default_obs = (0, 0, 0, 0)
        default_line = '\t'.join(str(x) for x in default_obs)
        # define the observation iterator
        if self.first == 'ignore':
            obs_it = self.gen_named_observations_unfilled(fin)
        else:
            obs_it = self.gen_named_observations_filled(fin)
        # yield chromosome names and observation lines
        for name, obs in obs_it:
            if obs is None:
                line = default_line
            else:
                line = '\t'.join(str(x) for x in obs)
            yield (name, line)


def get_form():
    """
    @return: the body of a form
    """
    sample_data = '\n'.join(g_sample_lines)
    form_objects = [
            Form.MultiLine('data_in', 'filtered pileup file', sample_data),
            Form.RadioGroup('first_info', 'first output position', [
                Form.RadioItem('f_0', '0'),
                Form.RadioItem('f_1', '1', True),
                Form.RadioItem('f_min', 'min'),
                Form.RadioItem('f_ignore', 'ignore')]),
            Form.RadioGroup('last_info', 'last output position', [
                Form.RadioItem('l_max', 'max', True),
                Form.RadioItem('l_1000', '1000'),
                Form.RadioItem('l_ignore', 'ignore')])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # unpack the first and last requested positions
    d_first = {'f_0':0, 'f_1':1, 'f_min':'min', 'f_ignore':'ignore'}
    d_last = {'l_max':'max', 'l_1000':1000, 'l_ignore':'ignore'}
    first = d_first[fs.first_info]
    last = d_last[fs.last_info]
    # create the scanner object which will be used for two passes
    scanner = Scanner(first, last, False)
    # Do the first pass; check for errors and gather chromosome info.
    names = set()
    fin = StringIO.StringIO(fs.data_in)
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
    out = StringIO.StringIO()
    print >> out, 'writing the first of', len(names), 'chromosomes:'
    print >> out
    fin = StringIO.StringIO(fs.data_in)
    for name, line in scanner.gen_named_lines(fin):
        if name == names[0]:
            print >> out, line
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
    if values[2] not in ambignt.g_resolve_nt:
        msg = 'the reference allele should be a nucleotide code: ' + values[2]
        raise Exception(msg)
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

def convert_row(row):
    """
    Return a flat tuple consisting of the reference and non-reference counts.
    The reference allele count is the first element of the tuple.
    The remaining allele counts are sorted in decreasing count order.
    @param row: a sequence of values in an expected format
    @return: a sufficient statistic
    """
    name, pos, ref = row[:3]
    A, C, G, T = row[6], row[8], row[10], row[12]
    acgt_counts = (A, C, G, T)
    nt_to_count = dict(zip('ACGT', acgt_counts))
    # hack the reference allele if it is ambiguous
    if ref not in list('ACGT'):
        nts = ambignt.g_resolve_nt[ref]
        count_nt_pairs = [(nt_to_count[nt], nt) for nt in nts]
        ref_count, ref = max(count_nt_pairs)
    # get the count of the reference allele followed by decreasing counts
    R = nt_to_count[ref]
    non_ref_counts = [nt_to_count[c] for c in 'ACGT' if c != ref]
    obs = [R] + list(reversed(sorted(non_ref_counts)))
    return tuple(obs)

def main(args):
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
    # create the scanner object which will be used for two passes
    scanner = Scanner(args.first, args.last)
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

def first_position(value):
    try:
        v = int(value)
    except ValueError, e:
        v = None
    if v is None:
        if value in ('min', 'drosophila', 'ignore'):
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
        if value in ('max', 'drosophila', 'truncated_drosophila', 'ignore'):
            return value
        else:
            raise TypeError()
    else:
        if v < 0:
            raise TypeError()
        else:
            return v

def check_position_requests(first, last):
    msg = "either both or neither of {first, last} should be 'ignore'"
    if sum(1 for x in (first, last) if x == 'ignore') == 1:
        raise ValueError(msg)
    msg_a = "if the last position is 'drosophila' or 'truncated_drosophila' "
    msg_b = "then the first position should be 'drosophila'"
    if last in ('drosophila', 'truncated_drosophila'):
        if first != 'drosophila':
            raise ValueError(msg)
    msg = 'the high position should not be lower than the low position'
    if first == 'drosophila' or type(first) == int:
        if type(last) == int:
            low = first
            if first == 'drosophila':
                low = 1
            high = last
            if high < low:
                raise ValueError(msg)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--profile', action='store_true',
            help='profile the script to look for slow spots'),
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--outdir', default=os.getcwd(),
            help='write the chromosome files to this directory')
    parser.add_argument('--out_prefix', default='chromosome.',
            help='prefix added to the chromosome name in the output filename')
    parser.add_argument('--out_suffix', default='.txt',
            help='suffix added to the chromosome name in the output filename')
    parser.add_argument('--first', default='drosophila', type=first_position,
            metavar='{<int>, min, drosophila, ignore}',
            help='the first position in a chromosome'),
    parser.add_argument('--last', default='drosophila', type=last_position,
            metavar='{<int>, max, drosophila, truncated_drosophila, ignore}',
            help='the last position in a chromosome'),
    args = parser.parse_args()
    if args.profile:
        profile.run('main(args)')
    else:
        main(args)
