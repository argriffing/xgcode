"""A wrapper for HyPhy.
"""

import unittest
import StringIO

import Util
import Newick


# This is example HyPhy output using a low verbosity level and custom fprintfs.
sample_hyphy_output = """
{
{             0.158}
{             0.326}
{             0.304}
{             0.212}
}

og Likelihood = -583.86147901356;done) LF Evals/Sec: 1693    CPU Load: 0.7533  
Shared Parameters:
eqFreqA=0.108328
eqFreqC=0.447118
eqFreqG=0.33229
eqFreqA2=0.216354
eqFreqC2=0.218499
eqFreqG2=0.272269
kappa=3.01447
P=0.456199
eqFreqT=1-eqFreqA-eqFreqC-eqFreqG=0.112264
eqFreqT2=1-eqFreqA2-eqFreqC2-eqFreqG2=0.292878

Tree givenTree=(((HUMAN:0.0759443,CHIMPANZEE:0.151889)Node2:0.607554,GORILLA:0.227833)Node1:0.53161,ORANGUTAN:0.303777,GIBBON:0.379722);
Tree otherTree=(((HUMAN:0.0930523,CHIMPANZEE:0.186105)Node2:0.744418,GORILLA:0.279157)Node1:0.651366,ORANGUTAN:0.372209,GIBBON:0.465261);
sf1=0.379722
sf2=0.465261

Check messages.log details of this run.
"""


class HyphyNamespace:

    def __init__(self):
        self.processed_lines = []

    def get_processed_lines(self):
        """
        This member function is for debugging.
        @return: a list of processed lines
        """
        return self.processed_lines

    def process_line(self, line):
        """
        Possibly add an attribute.
        The attribute can be a floating point value or a newick tree.
        """
        self.processed_lines.append(line)
        arr = [x.strip() for x in line.split('=')]
        if len(arr) > 1:
            var, value = arr[0], arr[-1]
            var_arr = var.split()
            # read the log likelihood
            if var == 'Log Likelihood' and value.endswith(';'):
                try:
                    fvalue = float(value[:-1])
                    setattr(self, 'lnL', fvalue)
                    return
                except ValueError:
                    pass
            # read a tree assignment
            if len(var_arr) == 2 and var_arr[0] == 'Tree':
                treename = var_arr[1]
                tree_string = value
                tree = Newick.parse(tree_string, Newick.NewickTree)
                setattr(self, treename, tree)
                return
            # read a standard variable assignment
            if len(var_arr) == 1:
                try:
                    fvalue = float(value)
                    setattr(self, var, fvalue)
                    return
                except ValueError:
                    pass


def get_hyphy_namespace(lines):
    """
    @param lines: lines of HyPhy output
    @return: a HyphyNamespace object
    """
    # process each line of the hyphy output
    ns = HyphyNamespace()
    for line in Util.stripped_lines(lines):
        ns.process_line(line)
    return ns


class TestHyphy(unittest.TestCase):

    def test_foo(self):
        get_hyphy_namespace(StringIO.StringIO(sample_hyphy_output))


def main():
    print '--test'

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestHyphy)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

