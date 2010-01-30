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


def gen_rows(line_source):
    for line in line_source:
        line = line.strip()
        if line:
            yield line.split()

def main(args):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile')
    parser.add_argument('--outdir')
    parser.add_argument('--force')
    args = parser.parse_args()
    main(args)
