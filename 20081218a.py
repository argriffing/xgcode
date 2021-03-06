"""
Get an amino acid Fasta alignment from some comma separated data.

The default tree is from the mm9_multiple_alignment page at UCSC.
"""

from StringIO import StringIO
import csv
import subprocess
import os

from SnippetUtil import HandlingError
import FelTree
import NewickIO
import Codon
import Fasta
import Util
import iterutils
import Form
import FormOut
import const

# TODO de-hardcode this path to the jar file
g_mapp_jar_path = '/home/argriffi/mapp-analysis/MAPP.jar'

g_default_csv_data = const.read('20100730t')

g_default_tree_string = const.read('20100730s')

def get_form():
    """
    @return: a list of form objects
    """
    default_tree_string = '(a:2, b:2, (c:2, d:2):2);'
    form_objects = [
            Form.MultiLine('data', 'raw data from a csv file',
                g_default_csv_data),
            Form.MultiLine('tree', 'a tree with branch lengths',
                g_default_tree_string)]
    return form_objects

def get_form_out():
    return FormOut.Fasta()

def get_amino_acid_alignment(table):
    """
    @param table: a table of data in some random format sent by Ferran Casals
    @return: a Fasta amino acid alignment object
    """
    if len(table) < 2:
        raise HandlingError('the data table should have at least two rows')
    first_row = table[0]
    if len(first_row) < 6:
        raise HandlingError(
                'the first row of the table has %d columns '
                'but at least six were expected' % len(first_row))
    if first_row[0].upper() != 'variant'.upper():
        raise HandlingError('expected the first column to be the variant')
    if first_row[1].upper() != 'chr'.upper():
        raise HandlingError('expected the second column to be the chromosome')
    if first_row[2].upper() != 'position'.upper():
        raise HandlingError('expected the third column to be the position')
    if first_row[3].upper() != 'Amino Acid Change'.upper():
        raise HandlingError(
                'expected the fourth column to be the amino acid change')
    if first_row[4].upper() != 'alleles'.upper():
        raise HandlingError(
                'expected the fifth column to be the nucleotide change')
    remaining_rows = table[1:]
    for row in remaining_rows:
        if len(row) != len(first_row):
            raise HandlingError(
                    'each row should have the same number of columns')
    # get the ordered taxa
    taxa = first_row[5:]
    if len(set(taxa)) != len(taxa):
        raise HandlingError('the same taxon appears in more than one column')
    # get the sequence of codons for each taxon
    codon_sequences = zip(*remaining_rows)[5:]
    # convert codon sequences to amino acid sequences
    aa_sequences = []
    for codon_sequence in codon_sequences:
        aa_list = []
        for codon in codon_sequence:
            codon = codon.upper()
            if codon == 'ND':
                aa = '-'
            elif codon in Codon.g_non_stop_codons:
                aa = Codon.g_codon_to_aa_letter[codon]
            elif codon in Codon.g_stop_codons:
                raise HandlingError(
                        'one of the codons is a stop codon: %s' % codon)
            else:
                raise HandlingError(
                        'one of the codons is invalid: %s' % codon)
            aa_list.append(aa)
        aa_sequences.append(''.join(aa_list))
    # return the alignment
    return Fasta.create_alignment(taxa, aa_sequences)

def get_alignment(data_string, tree_string):
    # convert the comma separated data into a table
    table = []
    for line in iterutils.stripped_lines(StringIO(data_string)):
        row = list(csv.reader(
            StringIO(line), delimiter=',', quotechar='"'))[0]
        table.append(row)
    # create the amino acid fasta alignment
    alignment = get_amino_acid_alignment(table)
    # create the tree
    tree = NewickIO.parse(tree_string, FelTree.NewickTree)
    # Make sure that the newick tree has all of the taxa
    # required by the alignment.
    tree_taxa_set = set(node.get_name() for node in tree.gen_tips())
    alignment_taxa_set = set(alignment.headers)
    weird_alignment_taxa = alignment_taxa_set - tree_taxa_set
    if weird_alignment_taxa:
        raise HandlingError(
                'the following taxa were not found '
                'in the tree: %s' % str(weird_taxa))
    # return the alignment
    return alignment

def get_mapp_output(data_string, tree_string):
    """
    @param data_string: a multi-line csv string
    @param tree_string: a multi-line newick string
    """
    # get the amino acid alignment
    alignment = get_alignment(data_string, tree_string)
    fasta_string = alignment.to_fasta_string()
    # make temporary files
    temp_fasta_filename = Util.create_tmp_file(fasta_string, suffix='.fa')
    temp_newick_filename = Util.create_tmp_file(tree_string, suffix='.tree')
    # call the mapp program
    cmd = [
            'gij',
            '-jar',
            g_mapp_jar_path,
            '-f',
            temp_fasta_filename,
            '-t',
            temp_newick_filename]
    # TODO use call or use communicate instead of wait?
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    process_output = p.stdout.read()
    process_error = p.stderr.read()
    # read the output
    return process_output

def get_response_content(fs):
    out = StringIO()
    # get the aligment
    alignment = get_alignment(fs.data, fs.tree)
    #print >> out, 'alignment:'
    print >> out, alignment.to_fasta_string()
    #print >> out
    # get the mapp output
    #mapp_output = get_mapp_output(fs.data, fs.tree, use_shell=True)
    #print >> out, 'mapp output:'
    #print >> out, mapp_output
    #print >> out
    # return the response
    """
    if fs.inline:
        response_headers = [('Content-Type', 'text/plain')]
        output_filename = 'raw-MAPP-output.txt'
    elif fs.attachment:
        response_headers = [('Content-Type', 'text/tab-separated-values')]
        output_filename = 'raw-MAPP-output.xls'
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.delivery, output_filename)))
    """
    return out.getvalue()
