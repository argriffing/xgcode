"""Given a fasta alignment, check it for syntax errors and briefly summarize it.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Fasta
import Form
import FormOut

def get_form():
    """
    @return: the body of a form
    """
    alignment_string = Fasta.brown_example_alignment.strip()
    return [Form.MultiLine('fasta', 'fasta alignment', alignment_string)]

def get_form_out():
    return FormOut.Report()

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    out = StringIO()
    try:
        alignment = Fasta.Alignment(StringIO(fs.fasta))
        print >> out, 'This is a valid alignment.'
    except Fasta.AlignmentError, e:
        alignment = None
        print >> out, 'This is not a valid alignment:', e
    if alignment:
        try:
            old_column_count = len(alignment.columns)
            alignment.force_nucleotide()
            removed_column_count = old_column_count - len(alignment.columns)
            if removed_column_count:
                print >> out, ('After removing %d' % removed_column_count),
                print >> out, 'columns this is a valid nucleotide alignment.'
            else:
                print >> out, 'This is a valid nucleotide alignment.'
        except Fasta.AlignmentError, e:
            print >> out, 'This is not a valid nucleotide alignment:', e
    for header, seq in Fasta.gen_header_sequence_pairs(StringIO(fs.fasta)):
        print >> out, '%s: %d' % (header, len(seq))
    response_headers = [('Content-Type', 'text/plain')]
    return response_headers, out.getvalue().strip()
