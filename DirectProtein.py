"""This module is related to the Direct Protein model.

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
import iterutils

# for testing
import PhyLikelihood
import Newick
import Fasta

from Codon import g_sorted_nt_letters as nt_ordered
from Codon import g_sorted_aa_letters as aa_ordered
from Codon import g_sorted_non_stop_codons as codons_ordered


# This tree is somewhat realistic but may not be properly scaled.
sample_tree_string = """
(
    ((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7,
    Orangutan:0.4,
    Gibbon:0.5
);
"""

# This codon alignment is unrealistic.
long_sample_codon_alignment_string = """
>Gibbon
TGTGAAGTGCTGCCTTTATTCTACTCCCCTTACCCTAGATTCATTCGTTTAGAGGGTGCC
CAGATCAGACGGTTCATAATAGGCCCGGAAGGCCACGGAGGACGATATCCACCCCTCTTC
CTATGCAATTCATGTCCTAGGGTCTCTTAT
>Orangutan
TGTGAAGTGGTGCCTCTATTCTACACCCCTTATCGTCGATTCATTAGTCTAGAAGGTGCC
CAAATTAGACGGTTCATAATAGGCCCGCGAGGGTGTGCAGGTCGATTTCCACCCCTATAT
TCGTGCAATTTCTCCCCTAGGGTCTCCTAT
>Gorilla
TATGATGTACCGCCTTTACTCTACATCCCACACGCGAAGCTCATTAGTCCAGAAGAGTCC
CAGAACAGAGTGGTCGTAAGGACCGCTGGAAGCCATGCAGGTAGATTCCCGCCCTTTACC
CCGTGCAATGTCCTCCCTATTGTTTCATAT
>Chimpanzee
TGCGAGAGGGCTCGTTTACTGTACATGCCTTGCTCTCAACCTACTAGTCCAGGGGGGCCC
TTTAACACATCTACCGTGATGAGCTCCGTATGGCATGCAGGTAGATTCCCGCTCTTAACC
TCTTGTAATGTCTTCCGTTTTGTTTCATAT
>Human
TGCGAGAGGGCTCCTTTACTGTACATGCCTTACTCTCAACCCATTAGTCCAGAGGGGCCC
CTTAGCAAACTTACCGTGATGAGCTCCGTTTGGCATGCAGGTAGATTCCCGCCCTTAACC
TCCTGTAATGTCTTCCGTTTCGTTTCATAT
"""

# This codon alignment is unrealistic.
short_sample_codon_alignment_string = """
>Gibbon
CCGTCGTACCGGCTGAACGGTTTCGATCGA
>Orangutan
TCGTCGAAACTGCTGAACGATGTTAATCAA
>Gorilla
TCCTCCCACCGGCGGAACGATATTAATAAA
>Chimpanzee
TCGTCCAGCCCGCGGAGGGCTATTCACAAA
>Human
TCGTCCAACCCGCGGAAGGCTATTCATAAA
"""


# for sanity checking only
eps = .000000001

def almost_equals(a, b):
    return abs(a-b) < eps

def get_nt_distribution_and_aa_energies(mutation_distribution, aa_distribution):
    """
    Use this function to guess the mutation distribution.
    If the output of this function is near the observed nucleotide distribution,
    then the input mutation distribution was almost correct.
    @param mutation_distribution: an ordered list of nucleotide frequencies defined by the mutation process
    @param aa_distribution: an ordered list of amino acid frequencies defined by both the mutation and selection process
    @return: the observed nucleotide distribution conditioned on the input distributions, and the amino acid energies
    """
    # do some error checking
    eps = 0.000000001
    if len(mutation_distribution) != 4:
        raise ValueError('expected four nucleotides')
    if len(aa_distribution) != 20:
        raise ValueError('expected twenty amino acids')
    if not almost_equals(sum(mutation_distribution), 1.0):
        raise ValueError('found a mutation distribution that does not sum to 1.0')
    if not almost_equals(sum(aa_distribution), 1.0):
        raise ValueError('found an amino acid distribution that does not sum to 1.0')
    for value in mutation_distribution:
        if almost_equals(value, 0):
            raise ValueError('each nucleotide should have a positive weight')
    for value in aa_distribution:
        if almost_equals(value, 0):
            raise ValueError('each amino acid should have a positive weight')
    # first get codon weights that depend only on the mutation distribution
    # and get the amino acid weights that depend only on the mutation distribution
    nt_to_weight = dict(zip(nt_ordered, mutation_distribution))
    aa_to_weight = dict(zip(aa_ordered, [0] * 20))
    codon_to_weight = {}
    for codon in codons_ordered:
        aa = Codon.g_codon_to_aa_letter[codon]
        weight = iterutils.product(nt_to_weight[nt] for nt in codon)
        codon_to_weight[codon] = weight
        aa_to_weight[aa] += weight
    # rescale codon and amino acid weights to sum to one
    total_weight = sum(aa_to_weight.values())
    for codon in codons_ordered:
        codon_to_weight[codon] /= total_weight
    for aa in aa_ordered:
        aa_to_weight[aa] /= total_weight
    # now find the amino acid exponentiated negative energies that scale the codons to the correct stationary distribution
    aa_to_exp_neg_energy = {}
    for aa, target_proportion in zip(aa_ordered, aa_distribution):
        aa_to_exp_neg_energy[aa] = target_proportion / aa_to_weight[aa]
    # now recalculate the codon weights to match the target aa distribution
    for codon in codons_ordered:
        aa = Codon.g_codon_to_aa_letter[codon]
        codon_to_weight[codon] *= aa_to_exp_neg_energy[aa]
    if not almost_equals(sum(codon_to_weight.values()), 1.0):
        raise HandlingError('final codon weights do not sum to 1.0')
    # get the final nucleotide weights
    nt_to_final_weight = dict(zip(nt_ordered, [0]*4))
    for codon, weight in codon_to_weight.items():
        for nt in codon:
            nt_to_final_weight[nt] += weight
    total_nt_weight = float(sum(nt_to_final_weight.values()))
    for nt in nt_ordered:
        nt_to_final_weight[nt] /= total_nt_weight
    if not almost_equals(sum(nt_to_final_weight.values()), 1.0):
        raise HandlingError('final nucleotide weights do not sum to 1.0')
    # get the final nucleotide list
    final_nucleotide_list = [nt_to_final_weight[nt] for nt in nt_ordered]
    # get the final amino acid list
    final_amino_acid_list = [-math.log(aa_to_exp_neg_energy[aa]) for aa in aa_ordered]
    # ok center the final amino acid list
    mean_energy = sum(final_amino_acid_list) / float(len(final_amino_acid_list))
    final_amino_acid_list = [energy - mean_energy for energy in final_amino_acid_list]
    # return the lists
    return final_nucleotide_list, final_amino_acid_list


class DirectProteinRateMatrix(RateMatrix.RateMatrix):
    """
    A Direct Protein codon substitution model.
    """

    def __init__(self, kappa, nucleotide_weights, amino_acid_energies):
        """
        @param kappa: the transition to transversion ratio of the nucleotide mutation process
        @param nucleotide_weights: an array of the ordered stationary frequencies of the nucleotide mutation process
        @param amino_acid_energies: an array of the ordered amino acid effects on protein energy
        """
        # validate the input
        if kappa < 0:
            raise ValueError('kappa must not be negative')
        if len(nucleotide_weights) != len(nt_ordered):
            raise ValueError('a weight must be specified for each nucleotide')
        for weight in nucleotide_weights:
            if weight < 0:
                raise ValueError('no nucleotide weight should be negative')
        if len(amino_acid_energies) != len(aa_ordered):
            raise ValueError('each energy list should specify an energy for each amino acid')
        # get each off-diagonal element of the rate matrix in convenient dictionary form
        aa_to_energy = dict(zip(aa_ordered, amino_acid_energies))
        nt_to_weight = dict(zip(nt_ordered, nucleotide_weights))
        codon_rate_matrix = {}
        for ca in codons_ordered:
            for cb in codons_ordered:
                rate = 0
                # if the codons differ at a nucleotide then the rate is nonzero
                if Util.hamming_distance(ca, cb) == 1:
                    # start multiplying together some factors to define the rate
                    rate = 1
                    # multiply by the factor due to nucleotide stationary frequency differences
                    for a, b in zip(ca, cb):
                        if a != b:
                            rate *= nt_to_weight[b]
                    # multiply by the factor due to transition / transversion rate
                    for a, b in zip(ca, cb):
                        if a != b:
                            if a+b in ('AG', 'GA', 'CT', 'TC'):
                                rate *= kappa
                    # multiply by the factor due to the amino acid energy difference
                    ea = aa_to_energy[Codon.g_codon_to_aa_letter[ca]]
                    eb = aa_to_energy[Codon.g_codon_to_aa_letter[cb]]
                    if ea != eb:
                        energy_difference = eb - ea
                        numerator = -energy_difference
                        denominator = 1 - math.exp(energy_difference)
                        rate *= numerator
                        rate /= denominator
                codon_rate_matrix[(ca, cb)] = rate
        # fill each diagonal element of the rate matrix
        for codon in codons_ordered:
            rate_away = sum(codon_rate_matrix[(codon, cb)] for cb in codons_ordered)
            codon_rate_matrix[(codon, codon)] = -rate_away
        # get the codon stationary distribution
        codon_to_stat_weight = {}
        for codon in codons_ordered:
            energy = aa_to_energy[Codon.g_codon_to_aa_letter[codon]]
            weight = 1
            weight *= math.exp(-energy)
            for nt in codon:
                weight *= nt_to_weight[nt]
            codon_to_stat_weight[codon] = weight
        # call the base class constructor
        row_major_rate_matrix = MatrixUtil.dict_to_row_major(codon_rate_matrix, codons_ordered, codons_ordered)
        RateMatrix.RateMatrix.__init__(self, row_major_rate_matrix, codons_ordered)
        # use a custom stationary state distribution without doing eigendecomposition
        total_weight = sum(codon_to_stat_weight.values())
        self.stationary_distribution = [codon_to_stat_weight[codon] / total_weight for codon in codons_ordered]
        # save the amino acid energies
        self.amino_acid_energies = amino_acid_energies

    def get_selection(self, ancestral_amino_acid, mutant_amino_acid):
        """
        Get the selection value of the new mutation given a population of the ancestral allele.
        This value will be positive when the new mutation is associated with less free energy.
        """
        aa_to_energy = dict(zip(aa_ordered, self.amino_acid_energies))
        return aa_to_energy[ancestral_amino_acid] - aa_to_energy[mutant_amino_acid]

    def get_codon_distribution(self):
        """
        @return: a codon stationary distribution dictionary defined by the rate matrix
        """
        return dict(zip(codons_ordered, self.stationary_distribution))

    def get_aa_distribution(self):
        """
        @return: an amino acid stationary distribution dictionary defined by the rate matrix
        """
        codon_distribution = self.get_codon_distribution()
        return Codon.codon_distribution_to_aa_distribution(codon_distribution)

    def get_nt_distribution(self):
        """
        @return: a nucleotide stationary distribution dictionary defined by the rate matrix
        """
        codon_distribution = self.get_codon_distribution()
        return Codon.codon_distribution_to_nt_distribution(codon_distribution)



class DirectProteinMixture(SubModel.MixtureModel):
    """
    A mixture of L{DirectProtein} codon substitution models.
    """

    def __init__(self, kappa, nucleotide_weights, mixture_weights, amino_acid_energy_lists):
        """
        @param kappa: the transition to transversion ratio of the nucleotide mutation process
        @param nucleotide_weights: the stationary distribution of the nucleotide mutation process
        @param mixture_weights: these mixing parameters are part of the selection process
        @param amino_acid_energy_lists: these amino acid energy lists are part of the selection process
        """
        # validate the mixture specific parts of the input
        if len(mixture_weights) != len(amino_acid_energy_lists):
            raise ValueError('the number of mixture weights must be the same as the number of energy lists')
        for weight in mixture_weights:
            if weight < 0:
                raise ValueError('no mixture weight should be negative')
        # create the rate matrices
        rate_matrices = []
        for amino_acid_energies in amino_acid_energy_lists:
            rate_matrix = DirectProteinRateMatrix(kappa, nucleotide_weights, amino_acid_energies)
            rate_matrices.append(rate_matrix)
        # call the base class constructor
        total_mixture_weight = float(sum(mixture_weights))
        mixture_distribution = [weight / total_mixture_weight for weight in mixture_weights]
        SubModel.MixtureModel.__init__(self, mixture_distribution, rate_matrices)
        # save some parameters so that the mixture can be saved as xml
        self.kappa = kappa
        self.nucleotide_weights = nucleotide_weights
        self.mixture_weights = mixture_weights
        self.amino_acid_energy_lists = amino_acid_energy_lists

    def get_codon_stationary_distribution(self):
        """
        @return: the stationary distribution of codons in the mixture model
        """
        return self.get_stationary_distribution()

    def get_aa_stationary_distribution(self):
        """
        @return: the stationary distribution of amino acids in the mixture model
        """
        aa_to_weight = dict((aa, 0) for aa in aa_ordered)
        for codon, proportion in zip(codons_ordered, self.get_codon_stationary_distribution()):
            aa = Codon.g_codon_to_aa_letter[codon]
            aa_to_weight[aa] += proportion
        total_weight = sum(aa_to_weight.values())
        return [aa_to_weight[aa] / total_weight for aa in aa_ordered]

    def get_nt_stationary_distribution(self):
        """
        @return: the stationary distribution of nucleotides in the mixture model
        """
        nt_to_weight = dict((nt, 0) for nt in nt_ordered)
        for codon, proportion in zip(codons_ordered, self.get_codon_stationary_distribution()):
            for nt in codon:
                nt_to_weight[nt] += proportion
        total_weight = sum(nt_to_weight.values())
        return [nt_to_weight[nt] / total_weight for nt in nt_ordered]

    def to_element_tree(self):
        """
        @return: an xml ElementTree representing the codon substitution mixture model
        """
        root = ET.Element('model')
        mutation = ET.SubElement(root, 'mutation')
        mutation.set('kappa', str(self.kappa))
        distribution = ET.SubElement(mutation, 'distribution')
        for nt, weight in zip(nt_ordered, self.nucleotide_weights):
            node = ET.SubElement(distribution, 'nt')
            node.set('symbol', nt)
            node.set('weight', str(weight))
        selection = ET.SubElement(root, 'selection')
        for mixture_weight, energy_list in zip(self.mixture_weights, self.amino_acid_energy_lists):
            category = ET.SubElement(selection, 'category')
            category.set('weight', str(mixture_weight))
            for aa, energy in zip(aa_ordered, energy_list):
                node = ET.SubElement(category, 'aa')
                node.set('symbol', aa)
                node.set('energy', str(energy))
        return ET.ElementTree(root)


class TestDirectProtein(unittest.TestCase):

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

    def test_likelihood_calculation(self):
        # get a tree
        tree = Newick.parse(sample_tree_string, Newick.NewickTree)
        # get a model
        input_xml_string = get_sample_xml_string()
        model = deserialize_mixture_model(input_xml_string)
        # get an alignment
        alignment = Fasta.CodonAlignment(StringIO(long_sample_codon_alignment_string))
        # get the likelihood
        log_likelihood = PhyLikelihood.get_log_likelihood(tree, alignment, model)

    def test_shooting_A(self):
        """
        Test the function that shoots towards the stationary nucleotide distribution.
        """
        expected_nt_dist = [
                0.21764166937,
                0.237427275677,
                0.290740146845,
                0.254190908108
                ]
        expected_centered_amino_acid_energies = [0.1] * 20
        expected_centered_amino_acid_energies[-2] = -1.9
        mut_nt_dist = [0.25, 0.25, 0.25, 0.25]
        aa_dist = [
                0.0593568189192,
                0.0296784094596,
                0.0296784094596,
                0.0296784094596,
                0.0296784094596,
                0.0593568189192,
                0.0296784094596,
                0.0445176141894,
                0.0296784094596,
                0.0890352283788,
                0.0148392047298,
                0.0296784094596,
                0.0593568189192,
                0.0296784094596,
                0.0890352283788,
                0.0890352283788,
                0.0593568189192,
                0.0593568189192,
                0.109647716212,
                0.0296784094596
                ]
        obs_nt_dist, obs_aa_energies = get_nt_distribution_and_aa_energies(mut_nt_dist, aa_dist)
        self.assertEquals(len(obs_nt_dist), len(expected_nt_dist))
        self.assertEquals(len(obs_aa_energies), len(expected_centered_amino_acid_energies))
        for observed, expected in zip(obs_nt_dist, expected_nt_dist):
            self.assertAlmostEquals(observed, expected)
        for observed, expected in zip(obs_aa_energies, expected_centered_amino_acid_energies):
            self.assertAlmostEquals(observed, expected)

    def test_shooting_B(self):
        """
        Test the function that shoots towards the stationary nucleotide distribution.
        """
        expected_nt_dist = [
                0.702432005526,
                0.106000805663,
                0.105569068691,
                0.0859981201201
                ]
        expected_centered_amino_acid_energies = [0.1] * 20
        expected_centered_amino_acid_energies[-2] = -1.9
        mut_nt_dist = [0.7, 0.1, 0.1, 0.1]
        aa_dist = [
                0.0106000805663,
                0.00212001611326,
                0.0148401127928,
                0.0593604511712,
                0.00212001611326,
                0.0106000805663,
                0.0148401127928,
                0.0667805075676,
                0.415523158198,
                0.0190801450193,
                0.0074200563964,
                0.10388078955,
                0.0106000805663,
                0.0593604511712,
                0.0699605317375,
                0.0254401933591,
                0.074200563964,
                0.0106000805663,
                0.00783245899575,
                0.0148401127928
                ]
        obs_nt_dist, obs_aa_energies = get_nt_distribution_and_aa_energies(mut_nt_dist, aa_dist)
        self.assertEquals(len(obs_nt_dist), len(expected_nt_dist))
        self.assertEquals(len(obs_aa_energies), len(expected_centered_amino_acid_energies))
        for observed, expected in zip(obs_nt_dist, expected_nt_dist):
            self.assertAlmostEquals(observed, expected)
        for observed, expected in zip(obs_aa_energies, expected_centered_amino_acid_energies):
            self.assertAlmostEquals(observed, expected)



def deserialize_mixture_model(xml_string):
    """
    @param xml_string: the xml string representing the substitution model
    @return: a L{DirectProteinMixture} object
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
    nucleotide_weights = [float(nt_to_weight[nt]) for nt in nt_ordered]
    # get the selection parameters
    selection = root.find('selection')
    mixture_weights = []
    amino_acid_energy_lists = []
    for category in selection:
        mixture_weights.append(float(category.get('weight')))
        aa_to_energy = {}
        for element in category:
            aa_to_energy[element.get('symbol')] = float(element.get('energy'))
        energy_list = [float(aa_to_energy[aa]) for aa in aa_ordered]
        amino_acid_energy_lists.append(energy_list)
    # create the mixture model object
    return DirectProteinMixture(kappa, nucleotide_weights, mixture_weights, amino_acid_energy_lists)

def get_sample_xml_string():
    """
    @return: a multi line xml string representing a L{DirectProteinMixture}.
    """
    root = ET.Element('model')
    mutation = ET.SubElement(root, 'mutation')
    mutation.set('kappa', '2.0')
    distribution = ET.SubElement(mutation, 'distribution')
    for nt in nt_ordered:
        node = ET.SubElement(distribution, 'nt')
        node.set('symbol', nt)
        node.set('weight', '1.0')
    selection = ET.SubElement(root, 'selection')
    for i in range(3):
        category = ET.SubElement(selection, 'category')
        category.set('weight', '3.0')
        for aa in aa_ordered:
            node = ET.SubElement(category, 'aa')
            node.set('symbol', aa)
            node.set('energy', '2.0')
    # modify the contents so that the tree is shown as indented
    XmlUtil.indent(root)
    # get the string representing the tree
    tree = ET.ElementTree(root)
    out = StringIO()
    tree.write(out)
    return out.getvalue()

def demo_xml():
    # create the tree
    root = ET.Element("html")
    head = ET.SubElement(root, "head")
    title = ET.SubElement(head, "title")
    title.text = "Page Title"
    body = ET.SubElement(root, "body")
    body.set("bgcolor", "#ffffff")
    body.text = "Hello, World!"
    tree = ET.ElementTree(root)
    # show the tree
    out = StringIO()
    tree.write(out)
    print out.getvalue()

def demo_codon_xml():
    print get_sample_xml_string()

if __name__ == '__main__':
    unittest.main()
