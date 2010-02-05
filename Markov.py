"""
Do some stuff with discrete time Markov chains.
"""

'http://www.gutenberg.org/files/98/98.txt'

from StringIO import StringIO
import unittest
import math

import Codon
import Util

class TestRateMatrix(unittest.TestCase):
    def test_foo(self):
        pass
    def test_bar(self):
        pass

def codon_rate_matrix_to_html_string(rate_matrix):
    codons = list(sorted(Codon.g_non_stop_codons))
    arr = []
    first_row_arr = ['<tr><td></td>'] + ['<td>%s</td>' % codon for codon in codons] + ['</tr>']
    arr.append(''.join(first_row_arr))
    for ca in Codon.g_non_stop_codons:
        row_arr = ['<tr><td>%s</td>' % ca] + ['<td>%.3f</td>' % rate_matrix[(ca, cb)] for cb in codons] + ['</tr>']
        arr.append(''.join(row_arr))
    return '\n'.join(arr)


def main():
    import CodonFrequency
    print '<html>'
    print '<head><style type="text/css">td{font-size:x-small;}</style></head>'
    print '<body>'
    for cf in (CodonFrequency.codon_frequency_a, CodonFrequency.codon_frequency_b):
        distribution = dict((codon, cf.codon_to_non_stop_proportion(codon)) for codon in Codon.g_non_stop_codons)
        rate_matrix = get_gy94_rate_matrix(distribution, 2, .01)
        table_guts = codon_rate_matrix_to_html_string(rate_matrix)
        print '<table>' + table_guts + '</table><br/><br/>'
    print '</body>'
    print '</html>'

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestRateMatrix)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()


