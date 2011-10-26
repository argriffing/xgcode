"""Convert multiple filtered pileup files to observation files for HMM.

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
from contextlib import nested
from multiprocessing import Pool
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import DGRP
import ambignt
import iterfiller
import const

#FIXME this might be the wrong file
g_sample_data = const.read('20100730b')

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
                Form.CheckItem('reqambig', 'o.o.b. refs must be N', True),
                Form.CheckItem('fill', 'fill with default values', True)])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # unpack the first and last requested positions
    low = {'low_0':0, 'low_1':1, 'low_none':None}[fs.low_info]
    high = {'high_1000':1000, 'high_none':None}[fs.high_info]
    # patch the fs object with the new values
    fs.low = low
    fs.high = high
    # process the lines according to the options
    out = StringIO()
    fin = StringIO(fs.data_in)
    for i, line in enumerate(gen_output_lines(fs, fin)):
        if i+1 >= 2000:
            msg = 'too many lines of output: ' + str(i+1)
            raise HandlingError(msg)
        print >> out, line
    return out.getvalue()

def gen_output_lines(args, fin):
    """
    Yield observation lines, given a filtered pileup file open for reading.
    """
    # unpack some relevant arguments
    reqambig = args.reqambig
    fill = args.fill
    errlow, errhigh = args.errlow, args.errhigh
    low = 1 if args.low == 'drosophila' else args.low
    high = args.high
    # create some state maintained across input lines
    filler = None
    chrom_name = None
    name_to_drosophila_length = dict(DGRP.g_chromosome_length_pairs)
    # define the default line to write
    default_obs = (0, 0, 0, 0)
    # process the input file line by line
    for line in fin:
        srow = line.split()
        if not srow:
            continue
        row = DGRP.filtered_pileup_row_to_typed(srow)
        obs = DGRP.filtered_pileup_typed_to_obs(row)
        name, pos, ref = row[:3]
        if filler is None:
            # set the chromosome name
            chrom_name = name
            # if appropriate, update the high value using the chrom name
            if args.high == 'drosophila':
                high = name_to_drosophila_length.get(name, None)
                if high is None:
                    raise Exception('invalid fly chromosome: ' + name)
            else:
                high = args.high
            # define the filler generator object
            filler = iterfiller.FillerGenerator(low, high,
                    fill, errlow, errhigh, default_obs)
        # check the chromosome name for consistency
        if name != chrom_name:
            msg_a = 'conflicting chromosome names: '
            msg_b = ' '.join((name, chrom_name))
            raise Exception(msg_a + msg_b)
        # check for reference nucleotide weirdness
        if reqambig:
            if not filler.check_bounds(pos):
                if ref != 'N':
                    msg_a = 'expected out of bounds reference nucleotides '
                    msg_b = 'to be N but found %s ' % ref
                    msg_c = 'at position %d of chrom %s' % (pos, name)
                    raise Exception(msg_a + msg_b + msg_c)
        # process lines emitted by the filler
        for value in filler.fill(pos, obs):
            yield '\t'.join(str(x) for x in value)
    # process final lines emitted by the filler
    for value in filler.finish():
        yield '\t'.join(str(x) for x in value)


def drosophila_position(value):
    """
    This is an argparse compatible type.
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

def convert_file_glom(glom):
    if len(glom) != 3:
        raise ValueError('expected three args glommed together')
    convert_file(*glom)

def convert_file(args, fpath_in, fpath_out):
    if args.dryrun:
        with open(fpath_in) as fin:
            for line in gen_output_lines(args, fin):
                pass
    else:
        with nested(open(fpath_in), open(fpath_out, 'w')) as (fin, fout):
            for line in gen_output_lines(args, fin):
                fout.write(line + '\n')

def main(args):
    # Define the output directory.
    output_directory = os.path.abspath(os.path.expanduser(args.outdir))
    # Make sure that the output directory exists.
    if not os.path.isdir(output_directory):
        if args.force:
            os.makedirs(output_directory)
        else:
            msg_a = 'not an existing directory: ' 
            msg_b = output_directory
            raise Exception(msg_a + msg_b)
    # Define the input and output files.
    fpaths_in = []
    fpaths_out = []
    for name in args.infiles:
        # get the full path to the input file
        fpath_in = os.path.abspath(os.path.expanduser(name))
        fpaths_in.append(fpath_in)
        # break the input file name into a directory part and a name part
        head_in, tail_in = os.path.split(fpath_in)
        # get the full path to the corresponding output file
        if head_in == output_directory:
            msg_a = 'input and output directories are the same: '
            msg_b = head_in
            raise Exception(msg_a + msg_b)
        fpath_out = os.path.join(output_directory, tail_in)
        fpaths_out.append(fpath_out)
    # Make that we are not overwriting anything important.
    if not args.force:
        for fpath_out in fpaths_out:
            if os.path.isfile(fpath_out):
                msg_a = 'would overwrite the file: '
                msg_b = fpath_out
                raise Exception(msg_a + msg_b)
    # Process each file using one or more processes.
    if args.nprocesses < 2:
        for fpath_in, fpath_out in zip(fpaths_in, fpaths_out):
            convert_file(args, fpath_in, fpath_out)
    else:
        glommed = [(args, fi, fo) for fi, fo in zip(fpaths_in, fpaths_out)]
        pool = Pool(processes=args.nprocesses)
        pool.map(convert_file_glom, glommed)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--nprocesses', default=1, type=int,
            help='use multiprocessing with a pool of this many processes')
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--dryrun', action='store_true',
            help='do not write any files')
    parser.add_argument('--reqambig', action='store_true',
            help='require ambig reference nt for out of bounds positions')
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
    parser.add_argument('infiles', nargs='+')
    args = parser.parse_args()
    main(args)
