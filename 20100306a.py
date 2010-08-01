"""Create a demo GFF track.
"""

from StringIO import StringIO

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.SingleLine('chromosome', 'chromosome name',
                '2L'),
            Form.Integer('npositions',
                'number of positions in the chromosome',
                23011546),
            Form.Integer('feature_width',
                'number of positions spanned by a feature',
                100000)]
    return form_objects

def get_form_out():
    return FormOut.Gff('demo')


def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    for i, line in enumerate(gen_cantor_track_lines(fs)):
        print >> out, line
        if i > 100:
            raise HandlingError('too many features')
    return [('Content-Type', 'text/plain')], out.getvalue().strip()

def change_base(n, base=3):
    while n:
        n, r = divmod(n, base)
        yield r

def gen_cantor_track_lines(args):
    yield 'track name=cantor description="test pattern"'
    nfeatures = args.npositions / args.feature_width
    for i in range(nfeatures):
        if 1 not in change_base(i):
            start = i * args.feature_width + 1
            stop = start + args.feature_width - 1
            arr = [
                    args.chromosome,
                    'myprogram', 'myfeature',
                    start, stop,
                    '.', '.', '.',
                    'mygroup']
            line = '\t'.join(str(x) for x in arr)
            yield line

def main(args):
    for line in gen_cantor_track_lines(args):
        print line

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--chromosome', default='2L',
            help='chromosome name')
    parser.add_argument('--npositions', default=23011546, type=int,
            help='number of positions in the chromosome')
    parser.add_argument('--feature_width', default=10000, type=int,
            help='number of positions spanned by a feature')
    args = parser.parse_args()
    main(args)
