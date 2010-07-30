"""Get an amino acid Fasta alignment from some comma separated data.

The default tree is from the mm9_multiple_alignment page at UCSC.
"""

from StringIO import StringIO
import csv
import tempfile
import subprocess
import os

from SnippetUtil import HandlingError
import FelTree
import NewickIO
import Codon
import Fasta
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
    """
            Form.RadioGroup('delivery', 'delivery', [
                Form.RadioItem('inline', 'view as text', True),
                Form.RadioItem('attachment', 'download as a tab separated data file')])]
    """
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
        msg_a = 'the first row of the table has %d columns ' % len(first_row)
        msg_b = 'but at least six were expected'
        raise HandlingError(msg_a + msg_b)
    if first_row[0].upper() != 'variant'.upper():
        raise HandlingError('expected the first column to be the variant')
    if first_row[1].upper() != 'chr'.upper():
        raise HandlingError('expected the second column to be the chromosome')
    if first_row[2].upper() != 'position'.upper():
        raise HandlingError('expected the third column to be the position')
    if first_row[3].upper() != 'Amino Acid Change'.upper():
        msg = 'expected the fourth column to be the amino acid change'
        raise HandlingError(msg)
    if first_row[4].upper() != 'alleles'.upper():
        msg = 'expected the fifth column to be the nucleotide change'
        raise HandlingError(msg)
    remaining_rows = table[1:]
    for row in remaining_rows:
        if len(row) != len(first_row):
            msg = 'each row should have the same number of columns'
            raise HandlingError(msg)
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
                msg = 'one of the codons is a stop codon: %s' % codon
                raise HandlingError(msg)
            else:
                msg = 'one of the codons is invalid: %s' % codon
                raise HandlingError(msg)
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
        msg_a = 'the following taxa were not found '
        msg_b = 'in the tree: %s' % str(weird_taxa)
        raise HandlingError(msg_a + msg_b)
    # return the alignment
    return alignment

def get_mapp_output(data_string, tree_string):
    """
    @param data_string: a multi-line csv string
    @param tree_string: a multi-line newick string
    """
    # get the amino acid alignment
    alignment = get_alignment(data_string, tree_string)
    # Get some temporary filenames for the alignment,
    # the tree, and the MAPP output.
    temp_fasta_filename = tempfile.mktemp(suffix='.fa')
    temp_newick_filename = tempfile.mktemp(suffix='.tree')
    # write the temporary fasta file
    temp_fasta_file = open(temp_fasta_filename, 'w')
    print >> temp_fasta_file, alignment.to_fasta_string()
    temp_fasta_file.close()
    # write the temporary newick tree file
    temp_newick_file = open(temp_newick_filename, 'w')
    print >> temp_newick_file, tree_string
    temp_newick_file.close()
    # call the mapp program
    cmd = [
            'gij',
            '-jar',
            g_mapp_jar_path,
            '-f',
            temp_fasta_filename,
            '-t',
            temp_newick_filename]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    process_output = p.stdout.read()
    process_error = p.stderr.read()
    # read the output
    return process_output


def get_response(fs):
    """
    @param fs: a decorated FieldStorage object
    @return: a (response_headers, response_text) pair
    """
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
    # write the response
    response_headers = [('Content-Type', 'text/plain')]
    """
    if fs.inline:
        response_headers = [('Content-Type', 'text/plain')]
        output_filename = 'raw-MAPP-output.txt'
    elif fs.attachment:
        response_headers = [('Content-Type', 'text/tab-separated-values')]
        output_filename = 'raw-MAPP-output.xls'
    response_headers.append(('Content-Disposition', "%s; filename=%s" % (fs.delivery, output_filename)))
    """
    return response_headers, out.getvalue().strip()

def main():
    print get_mapp_output(g_default_csv_data, g_default_tree_string)

if __name__ == '__main__':
    main()

