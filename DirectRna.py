"""This module is related to the Direct RNA model.

For more information about this model see:
"Population genetics without intraspecific data" by Thorne et al.
"""

from xml.etree import ElementTree as ET
import unittest
from StringIO import StringIO
import math

import Util
import Codon
import RateMatrix
import SubModel
import MatrixUtil
import XmlUtil


class DirectRnaRateMatrix(RateMatrix.RateMatrix):
    """
    A Direct RNA substitution model.
    """

    def __init__(self, kappa, nucleotide_weights, nucleotide_energies):
        """
        @param kappa: the transition to transversion ratio of the nucleotide mutation process
        @param nucleotide_weights: the stationary distribution of the nucleotide mutation process
        @param nucleotide_energies: selection presumably favors nucleotides associated with lower energies
        """
        # validate the input
        if kappa < 0:
            raise ValueError('kappa must not be negative')
        if len(nucleotide_weights) != len(Codon.g_nt_letters):
            raise ValueError('a weight must be specified for each nucleotide')
        for weight in nucleotide_weights:
            if weight < 0:
                raise ValueError('no nucleotide weight should be negative')
        if len(nucleotide_energies) != len(Codon.g_nt_letters):
            raise ValueError('each energy list should specify an energy for each nucleotide')
        # some utility variables
        ordered_nt_letters = list(sorted(Codon.g_nt_letters))
        # get each off-diagonal element of the rate matrix in convenient dictionary form
        nt_to_energy = dict(zip(ordered_nt_letters, nucleotide_energies))
        nt_to_weight = dict(zip(ordered_nt_letters, nucleotide_weights))
        nt_rate_matrix = {}
        for nta in ordered_nt_letters:
            for ntb in ordered_nt_letters:
                rate = 0
                if nta != ntb:
                    # consider the factor due to nucleotide stationary frequency differences
                    rate = nt_to_weight[ntb]
                    # multiply by the factor due to transition / transversion rate
                    if nta+ntb in ('AG', 'GA', 'CT', 'TC'):
                        rate *= kappa
                    # multiply by the factor due to the nucleotide energy difference
                    ea = nt_to_energy[nta]
                    eb = nt_to_energy[ntb]
                    if ea != eb:
                        energy_difference = eb - ea
                        numerator = -energy_difference
                        denominator = 1 - math.exp(energy_difference)
                        rate *= numerator
                        rate /= denominator
                nt_rate_matrix[(nta, ntb)] = rate
        # fill each diagonal element of the rate matrix
        for nt in ordered_nt_letters:
            rate_away = sum(nt_rate_matrix[(nt, ntb)] for ntb in ordered_nt_letters)
            nt_rate_matrix[(nt, nt)] = -rate_away
        # get the nucleotide stationary distribution
        nt_to_stat_weight = {}
        for nt in ordered_nt_letters:
            energy = nt_to_energy[nt]
            weight = 1
            weight *= math.exp(-energy)
            weight *= nt_to_weight[nt]
            nt_to_stat_weight[nt] = weight
        # call the base class constructor
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(nt_rate_matrix, ordered_nt_letters, ordered_nt_letters)
        RateMatrix.RateMatrix.__init__(self, row_major_rate_matrix, ordered_nt_letters)
        # use a custom stationary state distribution without doing eigendecomposition
        total_weight = sum(nt_to_stat_weight.values())
        self.stationary_distribution = [nt_to_stat_weight[nt] / total_weight for nt in ordered_nt_letters]
        # save the nucleotide energies
        self.nucleotide_energies = nucleotide_energies

    def get_selection(self, ancestral_nucleotide, mutant_nucleotide):
        """
        Get the selection value of the new mutation given a population of the ancestral allele.
        This value will be positive when the new mutation is associated with less free energy.
        """
        ordered_nt_letters = list(sorted(Codon.g_nt_letters))
        nt_to_energy = dict(zip(ordered_nt_letters, self.nucleotide_energies))
        nt_to_energy[ancestral_nucleotide] - aa_to_energy[mutant_nucleotide]


class DirectRnaMixture(SubModel.MixtureModel):
    """
    A mixture of L{DirectRna} nucleotide substitution models.
    """

    def __init__(self, kappa, nucleotide_weights, mixture_weights, nucleotide_energy_lists):
        """
        @param kappa: the transition to transversion ratio of the nucleotide mutation process
        @param nucleotide_weights: the stationary distribution of the nucleotide mutation process
        @param mixture_weights: these mixing parameters are part of the selection process
        @param nucleotide_energy_lists: these nucleotide energy lists are part of the selection process
        """
        # validate the mixture specific parts of the input
        if len(mixture_weights) != len(nucleotide_energy_lists):
            raise ValueError('the number of mixture weights must be the same as the number of energy lists')
        for weight in mixture_weights:
            if weight < 0:
                raise ValueError('no mixture weight should be negative')
        # create the rate matrices
        rate_matrices = []
        for nucleotide_energies in nucleotide_energy_lists:
            rate_matrix = DirectRnaRateMatrix(kappa, nucleotide_weights, nucleotide_energies)
            rate_matrices.append(rate_matrix)
        # call the base class constructor
        total_mixture_weight = float(sum(mixture_weights))
        mixture_distribution = [weight / total_mixture_weight for weight in mixture_weights]
        SubModel.MixtureModel.__init__(self, mixture_distribution, rate_matrices)
        # save some parameters so that the mixture can be saved as xml
        self.kappa = kappa
        self.nucleotide_weights = nucleotide_weights
        self.mixture_weights = mixture_weights
        self.nucleotide_energy_lists = nucleotide_energy_lists

    def get_nt_stationary_distribution(self):
        """
        @return: the stationary distribution of nucleotides in the mixture model
        """
        return self.get_stationary_distribution()

    def to_element_tree(self):
        """
        @return: an xml ElementTree representing the nucleotide substitution mixture model
        """
        root = ET.Element('model')
        mutation = ET.SubElement(root, 'mutation')
        mutation.set('kappa', str(self.kappa))
        distribution = ET.SubElement(mutation, 'distribution')
        for nt, weight in zip(Codon.g_nt_letters, self.nucleotide_weights):
            node = ET.SubElement(distribution, 'nt')
            node.set('symbol', nt)
            node.set('weight', str(weight))
        selection = ET.SubElement(root, 'selection')
        for mixture_weight, energy_list in zip(self.mixture_weights, self.nucleotide_energy_lists):
            category = ET.SubElement(selection, 'category')
            category.set('weight', str(mixture_weight))
            for nt, energy in zip(Codon.g_nt_letters, energy_list):
                node = ET.SubElement(category, 'nt')
                node.set('symbol', nt)
                node.set('energy', str(energy))
        return ET.ElementTree(root)


class TestDirectRna(unittest.TestCase):

    def test_sample_xml_string(self):
        """
        Verify that creating the sample xml string does not cause an exception.
        """
        # get the original sample xml string
        input_xml_string = get_sample_xml_string()
        # create a tree from the string
        element_tree = ET.parse(StringIO(input_xml_string))
        # create an xml string from the tree
        out = StringIO()
        element_tree.write(out)
        output_xml_string = out.getvalue()
        # verify that the output string is the same as the input string
        self.assertEquals(input_xml_string, output_xml_string)

    def test_serialization(self):
        """
        Verify that serialization and deserialization works.
        """
        # create the mixture model
        input_xml_string = get_sample_xml_string()
        mixture_model = deserialize_mixture_model(input_xml_string)
        # create an xml string from the mixture model
        element_tree = mixture_model.to_element_tree()
        XmlUtil.indent(element_tree.getroot())
        out = StringIO()
        element_tree.write(out)
        output_xml_string = out.getvalue()
        # verify that the xml string we get out is the same as the one we put in
        self.assertEquals(input_xml_string, output_xml_string)


def deserialize_mixture_model(xml_string):
    """
    @param xml_string: the xml string representing the substitution model
    @return: a L{DirectRnaMixture} object
    """
    element_tree = ET.parse(StringIO(xml_string))
    root = element_tree.getroot()
    # get the mutation parameters
    mutation = root.find('mutation')
    kappa = float(mutation.get('kappa'))
    distribution = mutation.find('distribution')
    nucleotide_weights = []
    nt_to_weight = {}
    for element in distribution:
        nt_to_weight[element.get('symbol')] = element.get('weight')
    nucleotide_weights = [float(nt_to_weight[nt]) for nt in Codon.g_nt_letters]
    # get the selection parameters
    selection = root.find('selection')
    mixture_weights = []
    nucleotide_energy_lists = []
    for category in selection:
        mixture_weights.append(float(category.get('weight')))
        nt_to_energy = {}
        for element in category:
            nt_to_energy[element.get('symbol')] = float(element.get('energy'))
        energy_list = [float(nt_to_energy[nt]) for nt in Codon.g_nt_letters]
        nucleotide_energy_lists.append(energy_list)
    # create the mixture model object
    return DirectRnaMixture(kappa, nucleotide_weights, mixture_weights, nucleotide_energy_lists)

def get_sample_xml_string():
    """
    @return: a multi line xml string representing a L{DirectRnaMixture}.
    """
    root = ET.Element('model')
    mutation = ET.SubElement(root, 'mutation')
    mutation.set('kappa', '2.0')
    distribution = ET.SubElement(mutation, 'distribution')
    for nt in Codon.g_nt_letters:
        node = ET.SubElement(distribution, 'nt')
        node.set('symbol', nt)
        node.set('weight', '1.0')
    selection = ET.SubElement(root, 'selection')
    for i in range(3):
        category = ET.SubElement(selection, 'category')
        category.set('weight', '3.0')
        for nt in Codon.g_nt_letters:
            node = ET.SubElement(category, 'nt')
            node.set('symbol', nt)
            node.set('energy', '2.0')
    # modify the contents so that the tree is shown as indented
    XmlUtil.indent(root)
    # get the string representing the tree
    tree = ET.ElementTree(root)
    out = StringIO()
    tree.write(out)
    return out.getvalue()

def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestDirectRna)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

