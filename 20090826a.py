"""Convert a resequencing file to a different format.

The input file has a comma separated header line like
'strain', 'chromosome', 'position', 'A', 'C', 'G', 'T', 'gap'
followed by conformant rows.
The output file will have tab separated data lines like
observed, a-count, c-count, g-count, t-count, gap-count
where 'observed' is 1 if data is available for the row
and is 0 if data is not available.
The position entries in the data rows of the input file are assumed to be
strictly monotonically increasing integers.
The number of rows in the output file will be
equal to (pmax - pmin + 1) where pmax is the largest position
entry in the input file, and pmin is the smallest position.
Rows in the output file that begin with 0 (that is, rows
corresponding to positions not represented in the input file)
will have undefined count entries.
They will probably be arbitrarily set to zero, but should be ignored.
The input file should represent a single chromosome and a single strain.
"""

import StringIO
import optparse
import sys

from SnippetUtil import HandlingError
import Form


class LineLimitError(Exception): pass


g_sample_input_rows = [
        ['strain', 'chromosome', 'position', 'A', 'C', 'G', 'T', 'gap'],
        ['222', '24', '101', '0', '0', '0', '0', '0'],
        ['222', '24', '103', '21', '1', '0', '0', '0'],
        ['222', '24', '108', '0', '12', '2', '0', '9'],
        ['222', '24', '109', '0', '0', '0', '0', '0'],
        ['222', '24', '110', '21', '1', '0', '0', '0'],
        ['222', '24', '122', '0', '12', '2', '666', '9']]


def get_form():
    """
    @return: the body of a form
    """
    sample_lines = [','.join(row) for row in g_sample_input_rows]
    form_objects = [
            Form.MultiLine('multiline_input', 'input file lines', '\n'.join(sample_lines))]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    line_source = StringIO.StringIO(fs.multiline_input)
    row_source = gen_rows(line_source)
    out = StringIO.StringIO()
    row_writer = RowWriter(out, line_limit=1000)
    try:
        process(row_source, row_writer)
    except LineLimitError, e:
        raise HandlingError(e)
    return [('Content-Type', 'text/plain')], out.getvalue().strip()

def gen_rows(line_source):
    for line in line_source:
        line = line.strip()
        if line:
            row = [x.strip() for x in line.split(',')]
            yield row


class RowWriter:

    def __init__(self, fout, line_limit=None):
        """
        @param fout: a file open for writing
        @param line_limit: limit the number of lines written
        """
        self.nrows_written = 0
        self.fout = fout
        self.line_limit = line_limit

    def write_row(self, row):
        """
        @param row: a sequence of strings
        """
        if self.line_limit is not None:
            if self.nrows_written == self.line_limit:
                raise LineLimitException('exceeded the limit of %d lines of output' % self.line_limit)
        print >> self.fout, '\t'.join(row)
        self.nrows_written += 1


def process(row_source, row_writer):
    """
    @param row_source: a source of data rows
    @param row_writer: write rows safely using this object
    """
    nlines_written = 0
    last_position = None
    header_row = next(row_source)
    for row in row_source:
        # read the row
        strain, chromosome, s_position, s_A, s_C, s_G, s_T, s_gap = row
        position = int(s_position)
        # write placeholder rows corresponding to skipped positions
        if last_position is not None:
            nskipped = position - last_position - 1
            for i in range(nskipped):
                row_writer.write_row(['0']*6)
        # update the last position
        last_position = position
        # write an informative row
        output_row = ['1', s_A, s_C, s_G, s_T, s_gap]
        row_writer.write_row(output_row)

def main(options, args):
    row_source = gen_rows(sys.stdin)
    row_writer = RowWriter(sys.stdout, line_limit=None)
    process(row_source, row_writer)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    options, args = parser.parse_args()
    main(options, args)
