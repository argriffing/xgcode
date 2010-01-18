"""Sample a nucleotide alignment given a tree, kappa, and a frequency mixture.

This is a silly model that I am using to experiment with HyPhy.
"""

import StringIO
import random

from xml.etree import ElementTree as ET

from SnippetUtil import HandlingError
import RateMatrix
import MatrixUtil
import XmlUtil
import SubModel
import Newick
import PhyLikelihood
import Nexus
import Form

def get_form():
    """
    @return: the body of a form
    """
    # define the tree string
    tree_string = '(((Human:0.1, Chimpanzee:0.2):0.8, Gorilla:0.3):0.7, Orangutan:0.4, Gibbon:0.5);'
    tree = Newick.parse(tree_string, Newick.NewickTree)
    formatted_tree_string = Newick.get_narrow_newick_string(tree, 60)
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'newick tree', formatted_tree_string),
            Form.MultiLine('model', 'DNA frequency mixture model', get_sample_xml_string().strip()),
            Form.Integer('ncols', 'sample this many nucleotide columns', 200, low=1, high=1000),
            Form.Integer('seed', 'random number seed', 314159, low=0),
            Form.RadioGroup('format', 'output format', [
                Form.RadioItem('fastaformat', 'fasta alignment'),
                Form.RadioItem('nexusformat', 'nexus alignment', True)])]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # parse the tree
    try:
        tree = Newick.parse(fs.tree, Newick.NewickTree)
        tree.assert_valid()
    except Newick.NewickSyntaxError, e:
        raise HandlingError(str(e))
    # get the normalized model
    mixture_model = deserialize_mixture_model(fs.model)
    # sample the alignment, possibly using a specified seed
    try:
        alignment = PhyLikelihood.simulate_alignment(tree, mixture_model, fs.ncols, fs.seed)
    except PhyLikelihood.SimulationError, e:
        raise HandlingError(e)
    # get the output string
    output_string = ''
    if fs.fastaformat:
        # the output is the alignment
        arr = []
        for node in tree.gen_tips():
            arr.append(alignment.get_fasta_sequence(node.name))
        alignment_string = '\n'.join(arr)
        output_string = alignment_string
    elif fs.nexusformat:
        # the output is the alignment and the tree
        nexus = Nexus.Nexus()
        nexus.tree = tree
        nexus.alignment = alignment
        output_string = str(nexus)
    # print the results
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, output_string

def get_sample_xml_string():
    """
    @return: a multi line xml string
    """
    # define the model
    kappa = 2
    category_weights = [1, 4, 5]
    nt_weights_list = [
            [1, 4, 4, 1],
            [2, 3, 3, 2],
            [1, 1, 1, 1]
            ]
    # create the model using the definition
    root = ET.Element('model')
    root.set('kappa', str(float(kappa)))
    for category_weight, nt_weights in zip(category_weights, nt_weights_list):
        category = ET.SubElement(root, 'category')
        category.set('weight', str(float(category_weight)))
        distribution = ET.SubElement(category, 'distribution')
        for nt, weight in zip(list('ACGT'), nt_weights):
            terminal = ET.SubElement(distribution, 'nt')
            terminal.set('symbol', nt)
            terminal.set('weight', str(float(weight)))
    # modify the contents so that the tree is shown as indented
    XmlUtil.indent(root)
    # get the string representing the tree
    tree = ET.ElementTree(root)
    out = StringIO.StringIO()
    tree.write(out)
    return out.getvalue()

def deserialize_mixture_model(xml_string):
    """
    Convert the xml string to a mixture model.
    @param xml_string: an xml string defining the mixture model
    @return: an unscaled mixture model object
    """
    # define the variables that define the model
    kappa = None
    category_weights = []
    nt_dicts = []
    # get the variables that define the model
    element_tree = ET.parse(StringIO.StringIO(xml_string))
    root = element_tree.getroot()
    kappa = float(root.get('kappa'))
    for category in root:
        category_weights.append(float(category.get('weight')))
        distribution = category.find('distribution')
        nt_dict = {}
        for terminal in distribution:
            nt_dict[terminal.get('symbol')] = float(terminal.get('weight'))
        total = sum(nt_dict.values())
        for nt in nt_dict:
            nt_dict[nt] /= total
        nt_dicts.append(nt_dict)
    # create a mixture model from the variables that define the model
    rate_matrix_objects = []
    for nt_dict in nt_dicts:
        rate_matrix_object = RateMatrix.get_unscaled_hky85_rate_matrix(nt_dict, kappa)
        rate_matrix_objects.append(rate_matrix_object)
    total = float(sum(category_weights))
    category_distribution = [weight / total for weight in category_weights]
    mixture_model = SubModel.MixtureModel(category_distribution, rate_matrix_objects)
    mixture_model.normalize()
    return mixture_model
