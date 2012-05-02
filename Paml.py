"""
A wrapper for the PAML tools of Ziheng Yang.
"""

import unittest
import subprocess
import popen2
import os

import Util
import Nexus
import Phylip
import Config
import iterutils

baseml_ctl = os.path.join(Config.data_path, 'baseml.ctl')
baseml_phylip = os.path.join(Config.data_path, 'baseml.phylip')
baseml_newick = os.path.join(Config.data_path, 'baseml.newick')
baseml_out = os.path.join(Config.data_path, 'baseml.out')

sample_paml_output = """
lnL(ntime:  0  np:  4):   -583.874748   +0.000000
  1.89613  0.21089  0.31930  0.17120

tree length =   3.00000

(((Human, Chimpanzee), Gorilla), Orangutan, Gibbon);

(((Human: 0.100000, Chimpanzee: 0.200000): 0.800000, Gorilla: 0.300000): 0.700000, Orangutan: 0.400000, Gibbon: 0.500000);

Detailed output identifying parameters
kappa under HKY85:  1.89613
base frequency parameters
  0.21089  0.31930  0.17120  0.29861
"""


def parse_hky_output(lines):
    """
    @param lines: lines of output
    @return: a dictionary with keys 'kappa', 'A', 'C', 'G', 'T', and 'lnL'
    """
    d = {}
    lines = Util.get_stripped_lines(lines)
    for line in lines:
        # read kappa
        if line.startswith('kappa under HKY85'):
            arr = [x.strip() for x in line.split(':')]
            d['kappa'] = float(arr[1])
        # read the log likelihood
        if line.startswith('lnL('):
            arr = line.split()
            d['lnL'] = float(arr[-2])
    # read the frequency parameters
    for first, second in iterutils.pairwise(lines):
        if first.startswith('base frequency parameters'):
            bases = list('TCAG')
            frequencies = [float(x) for x in second.split()]
            d.update(zip(bases, frequencies))
    return d


class PamlConfig:
    """
    Manages the configuration for the baseml PAML tool.
    """

    def __init__(self):
        self.config = {}

    def set_hky(self):
        """
        Set the configuration to the HKY model.
        """
        self.config = {
                # define the verbosity
                'noisy' : 0,
                'verbose' : 0,
                # define the input file names
                'seqfile' : baseml_phylip,
                'treefile' : baseml_newick,
                # define the output file name
                'outfile' : baseml_out,
                # use the HKY85 model
                'model' : 4,
                # unrooted
                'clock' : 0,
                # no topology inference
                'runmode' : 0,
                # no partitions
                'Mgene' : 0,
                # estimate kappa using an initial guess
                'fix_kappa' : 0,
                'kappa' : 1,
                # use a single rate for all sites
                'fix_alpha' : 1,
                'alpha' : 0,
                # use independent rates for sites
                'fix_rho' : 1,
                'rho' : 0,
                # estimate base frequencies
                'nhomo' : 1,
                # fix branch lengths
                'fix_blen' : 2
                }

    def to_ctl_string(self):
        arr = ['%s = %s' % item for item in self.config.items()]
        return '\n'.join(arr)



class TestPaml(unittest.TestCase):

    def test_foo(self):
        """
        Assert that creating a simple string does not fail.
        """
        config = PamlConfig()
        config.set_hky()
        config.to_ctl_string()


def run_hky(tree, alignment):
    """
    @param tree: a tree object
    @param alignment: an alignment object
    @return: messages from the program but not the results
    """
    # create the baseml.ctl control file
    config = PamlConfig()
    config.set_hky()
    config.to_ctl_string()
    with open(baseml_ctl, 'wt') as fout:
        print >> fout, config.to_ctl_string()
    # create the nexus object that defines the tree and alignment
    nexus = Nexus.get_sample_nexus_object()
    # create the baseml.newick tree file
    with open(baseml_newick, 'wt') as fout:
        print >> fout, nexus.tree.get_newick_string()
    # create the baseml.phylip alignment file
    phylip_string = Phylip.get_alignment_string_non_interleaved(nexus.alignment)
    with open(baseml_phylip, 'wt') as fout:
        print >> fout, phylip_string
    # change the current directory to the data directory
    with Util.remember_cwd():
        os.chdir(Config.data_path)
        # run PAML
        exe_path = Config.baseml_exe_path
        ctl_path = baseml_ctl
        #cmd = '%s %s > /dev/null' % (exe_path, ctl_path)
        #os.system(cmd)
        from_paml, to_paml = popen2.popen4([exe_path, ctl_path])
        #p = subprocess.Popen([cmd, arg], stdin=subprocess.PIPE, stdout=subprocess.PIPE, close_fds=True)
        # change back to the old directory
    return from_paml.read()

def main():
    nexus = Nexus.get_sample_nexus_object()
    run_hky(nexus.tree, nexus.alignment)
    with open(baseml_out) as fin:
        d = parse_hky_output(fin)
    print d

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestPaml)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

