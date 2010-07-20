"""Get the aligned amino acid region or column of a given genomic position.

Get the aligned amino acid region or column
associated with a given genomic position.
This snippet can also be run from the command line
to build the data directories required for its operation.
"""

from StringIO import StringIO
import sys
import os

from SnippetUtil import HandlingError
import KGEA
import Util
import Progress
from Form import RadioItem
import Form
import FormOut

# define some web locations that should probably be moved to a config file
g_base_dir = '/var/www/python_scripts/data/exon-alignments'
g_fasta_dir = g_base_dir + '/fasta-pieces'
g_index_dir = g_base_dir + '/piece-index-files'
g_valid_chromosome_strings_pathname = g_base_dir + '/valid-chromosome-strings.txt'


def get_form():
    """
    @return: the body of a form
    """
    # define the list of form objects
    form_objects = [
            Form.SingleLine('chromosome', 'chromosome', 'chr17'),
            Form.Integer('position', 'position', 70360012, low=0),
            Form.RadioGroup('output', 'output options', [
                RadioItem('show_alignment', 'show the aligned region', True),
                RadioItem('show_column', 'show the aligned column')])]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object decorated with field values
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    try:
        finder = KGEA.Finder(
                g_index_dir, g_valid_chromosome_strings_pathname, g_fasta_dir)
        if fs.show_alignment:
            lines = finder.get_alignment_lines(fs.chromosome, fs.position)
        elif fs.show_column:
            lines = finder.get_column_lines(fs.chromosome, fs.position)
        if lines:
            print >> out, '\n'.join(lines)
        else:
            msg = 'this position has no corresponding amino acid'
            raise KGEA.KGEAError(msg)
    except KGEA.KGEAError, e:
        raise HandlingError(e)
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue()


class MySyntaxError(Exception):
    pass


class MyConfigError(Exception):
    pass


def main(options, args):
    """
    @param options: from optparse
    @param args: from optparse
    @return: a response string
    """
    # get this file from the web directory
    # http://hgdownload.cse.ucsc.edu/goldenPath/hg18/multiz28way/alignments/
    original_fasta_filename = 'knownGene.exonAA.fa'
    # the original fasta file is broken into a bunch of pieces and put into this directory
    pieces_directory = 'fasta'
    # the index files that map genomic locations to fasta subfiles are in this directory
    index_directory = 'index'
    # this file keeps a list of valid chromosome names from the original fasta file
    chromosome_filename = 'chromosomes.txt'
    # assert that a command was given with the script
    if not args:
        raise MySyntaxError('no command was given')
    # try to dispatch the command
    command = args[0]
    command_args = args[1:]
    if command == 'split':
        if command_args:
            raise MySyntaxError('the split command does not take any arguments')
        # assert that the fasta directory has been created
        if not os.path.isdir(pieces_directory):
            err_lines = [
                    'The directory for the split fasta files was not found: ' + pieces_directory,
                    'Please create this directory or cd to its parent directory.']
            raise MyConfigError('\n'.join(err_lines))
        # assert that the current directory has the original huge fasta file
        pathnames = os.listdir('.')
        if original_fasta_filename not in pathnames:
            err_lines = [
                    'The file %s was not found in the current directory.' % original_fasta_filename,
                    'Please download this file from:',
                    'http://hgdownload.cse.ucsc.edu/goldenPath/hg18/multiz28way/alignments/']
            raise MyConfigError('\n'.join(err_lines))
        splitter = KGEA.Splitter(original_fasta_filename, pieces_directory)
        splitter.run(verbose=options.verbose)
        return ''
    elif command == 'index':
        if command_args:
            raise MySyntaxError('the index command does not take any arguments')
        # assert that the fasta directory has been created
        if not os.path.isdir(pieces_directory):
            err_lines = [
                    'The directory for the split fasta files was not found: ' + pieces_directory,
                    'If this directory exists somewhere else, then cd to its parent directory.',
                    'If this directory has not been created, then create it and run the split command.']
            raise MyConfigError('\n'.join(err_lines))
        # assert that the index directory has been created
        if not os.path.isdir(index_directory):
            err_lines = [
                    'The directory for the index files was not found: ' + index_directory,
                    'Please create this directory or cd to its parent directory.']
            raise MyConfigError('\n'.join(err_lines))
        indexer = KGEA.Indexer(index_directory, chromosome_filename, pieces_directory)
        indexer.run(verbose=options.verbose)
        return ''
    elif command == 'find-alignment':
        if len(command_args) != 2:
            raise MySyntaxError('the find-alignment command takes two arguments')
        # define the chromosome string and the chromosome position
        chromosome_string, chromosome_position_string = command_args
        # initialize the chromosome position and assert that it is plausible
        try:
            chromosome_position = int(chromosome_position_string)
        except ValueError, e:
            raise MySyntaxError('the chromosome position should be an integer')
        # assert that the fasta directory has been created
        if not os.path.isdir(pieces_directory):
            err_lines = [
                    'The directory for the split fasta files was not found: ' + pieces_directory,
                    'If this directory exists somewhere else, then cd to its parent directory.',
                    'If this directory has not been created, then create it and run the split command.']
            raise MyConfigError('\n'.join(err_lines))
        # assert that the index directory has been created
        if not os.path.isdir(index_directory):
            err_lines = [
                    'The directory for the index files was not found: ' + index_directory,
                    'If this directory exists somewhere else, then cd to its parent directory.',
                    'If this directory has not been created, then create it and run the index command.']
            raise MyConfigError('\n'.join(err_lines))
        # look for the alignment using the finder
        finder = KGEA.Finder(index_directory, chromosome_filename, pieces_directory)
        fasta_lines = finder.get_alignment_lines(chromosome_string, chromosome_position, verbose=options.verbose)
        if not fasta_lines:
            return 'no amino acid was found at this position'
        return '\n'.join(fasta_lines)
    elif command == 'find-column':
        if len(command_args) != 2:
            raise MySyntaxError('the find-column command takes two arguments')
        # define the chromosome string and the chromosome position
        chromosome_string, chromosome_position_string = command_args
        # initialize the chromosome position and assert that it is plausible
        try:
            chromosome_position = int(chromosome_position_string)
        except ValueError, e:
            raise MySyntaxError('the chromosome position should be an integer')
        # assert that the fasta directory has been created
        if not os.path.isdir(pieces_directory):
            err_lines = [
                    'The directory for the split fasta files was not found: ' + pieces_directory,
                    'If this directory exists somewhere else, then cd to its parent directory.',
                    'If this directory has not been created, then create it and run the split command.']
            raise MyConfigError('\n'.join(err_lines))
        # assert that the index directory has been created
        if not os.path.isdir(index_directory):
            err_lines = [
                    'The directory for the index files was not found: ' + index_directory,
                    'If this directory exists somewhere else, then cd to its parent directory.',
                    'If this directory has not been created, then create it and run the index command.']
            raise MyConfigError('\n'.join(err_lines))
        # look for the column using the finder
        finder = KGEA.Finder(index_directory, chromosome_filename, pieces_directory)
        column_lines = finder.get_column_lines(chromosome_string, chromosome_position, verbose=options.verbose)
        if not column_lines:
            return 'no amino acid was found at this position'
        return '\n'.join(column_lines)
    elif command == 'summarize':
        if command_args:
            raise MySyntaxError('the summarize command does not take any arguments')
        # assert that the current directory has the original huge fasta file
        pathnames = os.listdir('.')
        if original_fasta_filename not in pathnames:
            err_lines = [
                    'The file %s was not found in the current directory.' % original_fasta_filename,
                    'Please download this file from:',
                    'http://hgdownload.cse.ucsc.edu/goldenPath/hg18/multiz28way/alignments/']
            raise MyConfigError('\n'.join(err_lines))
        # initialize the progress bar
        nbytes_total = os.path.getsize(original_fasta_filename)
        pbar = Progress.Bar(nbytes_total)
        # initialize the summary
        mod3 = {0:0, 1:0, 2:0}
        length_diff_dict = {}
        # summarize by reading each alignment from the file
        approx_nbytes_read = 0
        fin = open(original_fasta_filename)
        for lines in Util.gen_paragraphs(fin):
            # process the lines
            header_line = lines[0]
            p = KGEA.LocationParser(header_line)
            genomic_length = (p.last_index - p.first_index) + 1
            mod3[genomic_length % 3] += 1
            diff = 3*p.length - genomic_length
            if diff not in length_diff_dict:
                length_diff_dict[diff] = 0
            length_diff_dict[diff] += 1
            # update the progress bar
            approx_nbytes_read += sum(len(line) for line in lines)
            pbar.update(approx_nbytes_read)
        fin.close()
        # finish the progress bar
        pbar.update(nbytes_total)
        # return the summary
        summary_lines = []
        summary_lines += ['genomic span of %d mod 3: %d sequences' % (i, mod3[i]) for i in range(3)]
        summary_lines.append('histogram of 3*aa_length - genomic span:')
        for key, value in sorted(length_diff_dict.items()):
            summary_lines.append('%d : %d' % (key, value))
        return '\n'.join(summary_lines)
    else:
        raise MySyntaxError('invalid command: ' + command)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    options, args = parser.parse_args()
    try:
        response_string = main(options, args)
        print response_string
    except MySyntaxError, e:
        print 'Syntax error.'
        print e
        print 'Here are some examples:'
        print 'python %s find-alignment chr7 15393678' % sys.argv[0]
        print 'python %s find-column chr7 15393678' % sys.argv[0]
        print 'python %s split' % sys.argv[0]
        print 'python %s index' % sys.argv[0]
        print 'python %s summarize' % sys.argv[0]
    except MyConfigError, e:
        print 'Configuration error.'
        print e
    except KGEA.KGEAError, e:
        print 'Error.'
        print e


