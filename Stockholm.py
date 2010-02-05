"""Write alignment file using the stockholm format.
"""

import unittest
from StringIO import StringIO

import Fasta
import Newick
import Monospace
import Util

# This example is from the RSmatch software documentation.
stockholm_sample_string = """
# STOCKHOLM 1.0
#=GF ID CBS
#=GF AC PF00571
#=GF DE CBS domain
#=GF AU Bateman A
#=GF CC CBS domains are small intracellular modules mostly found  
#=GF CC in 2 or four copies within a protein. 
#=GF SQ 67
#=GS O31698/18-71 AC O31698
#=GS O83071/192-246 AC O83071
#=GS O83071/259-312 AC O83071
#=GS O31698/88-139 AC O31698
#=GS O31698/88-139 OS Bacillus subtilis
O83071/192-246          MTCRAQLIAVPRASSLAE..AIACAQKM....RVSRVPVYERS
#=GR O83071/192-246 SA  999887756453524252..55152525....36463774777
O83071/259-312          MQHVSAPVFVFECTRLAY..VQHKLRAH....SRAVAIVLDEY
#=GR O83071/259-312 SS  CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEE
O31698/18-71            MIEADKVAHVQVGNNLEH..ALLVLTKT....GYTAIPVLDPS
#=GR O31698/18-71 SS    CCCHHHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEHHH
O31698/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31698/88-139 SS   CCCCCCCHHHHHHHHHHH..HEEEEEEE....EEEEEEEEEEH
#=GC SS_cons            CCCCCHHHHHHHHHHHHH..EEEEEEEE....EEEEEEEEEEH
O31699/88-139           EVMLTDIPRLHINDPIMK..GFGMVINN......GFVCVENDE
#=GR O31699/88-139 AS   ________________*__________________________
#=GR_O31699/88-139_IN   ____________1______________2__________0____
//
"""


class StockholmError(Exception):
    pass


class Stockholm:
    """
    As far as I am concerned, this is a data structure that aggregates an alignment and a tree.
    Set the alignment and the tree manually.
    """

    def __init__(self):
        # The stockholm object should have an alignment and a tree.
        self.alignment = None
        self.tree = None
        self.comments = []
        self.column_annotations = []

    def add_comment(self, comment):
        self.comments.append(comment)

    def add_column_annotation(self, annotation_name, annotation_value):
        self.column_annotations.append((annotation_name, annotation_value))

    def __str__(self):
        out = StringIO()
        alignment = self.alignment
        tree = self.tree
        # write the format identifier
        print >> out, '# STOCKHOLM 1.0'
        # write the comments
        for comment in self.comments:
            print >> out, '#=GF CC', comment
        # write the tree
        print >> out, '#=GF NH', tree.get_newick_string()
        # determine the max header length
        augmented_column_annotation_names = ['#=GC '+name for name, value in self.column_annotations]
        headers = alignment.headers + augmented_column_annotation_names
        max_header_length = max(len(header) for header in headers)
        # write the alignment
        for header, sequence in zip(alignment.headers, alignment.sequences):
            left_justified_header = Monospace.left_justify(header, max_header_length, ' ')
            print >> out, '%s %s' % (left_justified_header, sequence)
        # write the column annotations
        column_annotation_values = [value for name, value in self.column_annotations]
        for name, value in zip(augmented_column_annotation_names, column_annotation_values):
            justified_name = Monospace.left_justify(name, max_header_length, ' ')
            print >> out, justified_name, value
        # write the format terminator
        print >> out, '//'
        return out.getvalue()


class TestStockholm(unittest.TestCase):

    def test_stockholm(self):
        pass


def main():
    print '...'

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestStockholm)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()
