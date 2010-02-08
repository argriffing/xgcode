"""Get an amino acid Fasta alignment from some comma separated data.

The default tree is from the mm9_multiple_alignment page at UCSC.
"""

from StringIO import StringIO
import csv
import tempfile
import subprocess
import os

from SnippetUtil import HandlingError
import Config
import FelTree
import NewickIO
import Codon
import Fasta
import iterutils
import Form

# this is some sample data
g_default_csv_lines = [
        '"variant","chr","position","Amino Acid Change","alleles","Hg18","bosTau3","canFam2","danRer5","fr1","galGal3","gasAcu1","mm9","monDom4","ornAna1","oryLat1","panTro2","rheMac2","rn4","tetNig1"',
        '"v47291",17,70360012,"L-Q","T-A","CTG","TTG","CTG","CTG","TTA","CTG","CTG","TTG","CTG","ND","CTC","CTG","CTC","CTG","TTA"',
        '"v47321",17,70360024,"A-E","C-A","GCG","GCG","GCG","AGG","TCT","GGG","GCC","GCA","GAG","ND","TCC","GCG","GAG","GCA","GCC"',
        '"v47278",17,70362556,"H-R","A-G","CAC","CAT","CAC","CAC","CAC","CAC","CAT","CAC","CAT","ND","CAT","CGC","CAC","CAC","CAC"',
        '"v47269",17,70362745,"Q-H","G-T","CAG","CAG","CAG","CCC","CCC","CCG","CCC","CAG","CAG","ND","CCC","CAG","CAG","CAG","CCC"',
        '"v49364",7,15393678,"F-L","T-C","TTT","GTT","GTT","TTA","ND","GCT","ND","CTT","CAT","ND","ND","TTT","TTT","CTC","ND"']

# this is from genomewiki.ucsc.edu/index.php/Mm9_multiple_alignment
g_default_tree_lines = [
        '((((((((',
        '(((mm9:0.076274,rn4:0.084383):0.200607,',
        'cavPor2:0.202990):0.034350,',
        'oryCun1:0.208548):0.014587,',
        '((((((Hg18:0.005873,panTro2:0.007668):0.013037,',
        'ponAbe2:0.02):0.013037,rheMac2:0.031973):0.0365,',
        'calJac1:0.07):0.0365,otoGar1:0.151185):0.015682,',
        'tupBel1:0.162844):0.006272):0.019763,',
        '((sorAra1:0.248532,eriEur1:0.222255):0.045693,',
        '(((canFam2:0.101137,felCat3:0.098203):0.048213,',
        'equCab1:0.099323):0.007287,',
        'bosTau3:0.163945):0.012398):0.018928):0.030081,',
        '(dasNov1:0.133274,(loxAfr1:0.103030,',
        'echTel1:0.232706):0.049511):0.008424):0.213469,',
        'monDom4:0.320721):0.088647,',
        'ornAna1:0.488110):0.118797,',
        '(galGal3:0.395136,anoCar1:0.513962):0.093688):0.151358,',
        'xenTro2:0.778272):0.174596,',
        '(((tetNig1:0.203933,fr1:0.239587):0.203949,',
        '(gasAcu1:0.314162,oryLat1:0.501915):0.055354):0.346008,',
        'danRer5:0.730028):0.174596);']

def get_form():
    """
    @return: a list of form objects
    """
    default_tree_string = '(a:2, b:2, (c:2, d:2):2);'
    form_objects = [
            Form.MultiLine('data', 'raw data from a csv file',
                '\n'.join(g_default_csv_lines)),
            Form.MultiLine('tree', 'a tree with branch lengths',
                '\n'.join(g_default_tree_lines))]
    """
            Form.RadioGroup('delivery', 'delivery', [
                Form.RadioItem('inline', 'view as text', True),
                Form.RadioItem('attachment', 'download as a tab separated data file')])]
    """
    return form_objects

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
                msg = 'one of the codons is invalid: %s' % codon)
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
    mapp_jar_pathname = Config.mapp_exe_path + '/' + 'MAPP.jar'
    cmd = [
            'gij',
            '-jar',
            mapp_jar_pathname,
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
    data_string = '\n'.join(g_default_csv_lines)
    tree_string = '\n'.join(g_default_tree_lines)
    print get_mapp_output(data_string, tree_string)

if __name__ == '__main__':
    main()

