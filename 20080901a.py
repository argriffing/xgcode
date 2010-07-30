"""Summarize MAPP results.

MAPP stands for Multivariate Analysis of Protein Polymorphism.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import Codon
import Util
import Form
import FormOut
import const

g_mapp_output = const.read('20100730v')


class ColumnDataError(Exception):
    pass


class ColumnData:
    """
    Define data associated with the column of an alignment.
    """

    def __init__(self, gene, offset, wild, mutant):
        """
        @param gene: the name of a gene
        @param offset: the integer offset indexed starting at 1
        @param wild: the wild type amino acid
        @param mutant: the mutant type amino acid
        """
        if type(offset) is not type(1):
            raise ColumnDataError('invalid offset')
        if offset < 1:
            raise ColumnDataError('the offset must be greater or equal to one')
        if wild not in Codon.g_aa_letters:
            raise ColumnDataError('invalid wild type')
        if mutant not in Codon.g_aa_letters:
            raise ColumnDataError('invalid mutant type')
        self.gene = gene
        self.offset = offset
        self.wild = wild
        self.mutant = mutant

    def __str__(self):
        return self.gene + ':' + self.wild + str(self.offset) + self.mutant


def get_form():
    """
    @return: the body of a form
    """
    # define some default data
    column_info_list = (
            ColumnData('TMEM195', 279, 'F', 'L'),
            ColumnData('SHANK3', 552, 'R', 'W'),
            ColumnData('ARHGAP6', 76, 'G', 'D'),
            ColumnData('GRPR', 261, 'R', 'L'),
            ColumnData('DRP2', 203, 'T', 'M'),
            ColumnData('PLXNB3', 981, 'R', 'H'))
    data_lines = [str(info) for info in column_info_list]
    # define the form objects
    form_objects = [
            Form.MultiLine('mapp', 'MAPP output',
                g_mapp_output.strip()),
            Form.MultiLine('headers', 'alignment column headers',
                '\n'.join(data_lines))]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # read the headers
    headers = list(Util.stripped_lines(StringIO(fs.headers)))
    # read the tab separated MAPP output
    tsv_lists = []
    for line in StringIO(fs.mapp):
        if line.strip():
            tsv_list = [element.strip() for element in line.split('\t')]
            tsv_lists.append(tsv_list)
    # validate the headers, possibly providing default values
    if headers:
        if len(headers) != len(tsv_lists) - 1:
            msg_a = 'the number of headers should be one fewer than '
            msg_b = 'the number of MAPP lines'
            raise HandlingError(msg_a + msg_b)
    else:
        headers = [str(i+1) for i in range(len(tsv_lists) - 1)]
    # check input consistency
    length_set = set(len(tsv_list) for tsv_list in tsv_lists)
    if length_set != set([54]):
        msg_a = 'each line in the MAPP output should have 54 '
        msg_b = 'tab separated values: %s' % str(length_set)
        raise HandlingError(msg_a + msg_b)
    # read the p-values
    pvalue_lists = []
    for tsv_list in tsv_lists[1:]:
        pvalue_list = []
        for element in tsv_list[32:52]:
            try:
                pvalue = float(element)
            except ValueError:
                pvalue = None
            pvalue_list.append(pvalue)
        pvalue_lists.append(pvalue_list)
    # define the response
    out = StringIO()
    for header, pvalue_list in zip(headers, pvalue_lists):
        # if a p-value is bad then do not show anything for this column
        if None in pvalue_list:
            continue
        # show column information
        print >> out, header
        for pvalue, aa_letter in sorted(zip(pvalue_list, Codon.g_aa_letters)):
            print >> out, aa_letter, pvalue
        print >> out, ''
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
