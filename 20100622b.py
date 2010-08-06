"""Fungus tsv janitor.

Clean the fungus tsv file using ad hoc codes.
The output is an R table.
"""

from StringIO import StringIO

from SnippetUtil import HandlingError
import Form
import FormOut
import Util
import Carbone

g_tags = ['pca:convert']

g_info_rows = [
        ['IC', 'Location', 'Temp_C', 'Precip_mm', 'Species'],
        ['IC13', 'GA', '15', '600', 'Ap'],
        ['IC14', 'GA', '15', '600', 'Ap'],
        ['IC1493', 'LA', '15', '600', 'Ano'],
        ['IC1494', 'LA', '15', '600', 'Ano']]

g_info_lines = ['\t'.join(row) for row in g_info_rows]

g_input_headers = [
        'otu', 'location', 'temperature', 'precipitation', 'species']

g_output_headers = [
        'temperature', 'precipitation', 'otu']


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('info', 'a table with column but not row headers',
                '\n'.join(g_info_lines)),
            Form.MultiLine('input_headers', 'renamed column headers',
                '\n'.join(g_input_headers)),
            Form.MultiLine('output_headers', 'ordered output headers',
                '\n'.join(g_output_headers)),
            Form.RadioGroup('mdefgroup', 'missing input data', [
                Form.RadioItem('star_missing_in',
                    '* is interpreted as missing data', True),
                Form.RadioItem('NULL_missing_in',
                    'NULL is interpreted as missing data'),
                Form.RadioItem('no_missing_in',
                    'no data is missing')]),
            Form.RadioGroup('mhandlegroup', 'missing output data', [
                Form.RadioItem('remove_missing_out',
                    'remove rows that have missing data', True),
                Form.RadioItem('NA_missing_out',
                    'show missing data as NA')]),
            Form.CheckGroup('cleangroup', 'more options', [
                Form.CheckItem('add_indices',
                    'add row indices for R table compatibility', True),
                Form.CheckItem('clean_isolates',
                    'force first-column elements to be IC-prefixed', True)]),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.RTable('out')

def get_response_content(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    return process(fs, fs.info.splitlines(), fs.input_headers.splitlines(),
            fs.output_headers.splitlines()) + '\n'

def process(args, raw_info_lines, raw_input_headers, raw_output_headers):
    info_lines = Util.get_stripped_lines(raw_info_lines)
    rows = [line.split() for line in info_lines]
    # the number of columns should be consistent among rows
    if len(set(len(row) for row in rows)) != 1:
        msg = 'the number of columns should be consistent among rows'
        raise ValueError(msg)
    # break the list of rows into a header row and data rows
    header, data_rows = rows[0], rows[1:]
    # account for missing input data
    if args.star_missing_in:
        data_rows = [[None if v=='*' else v for v in r] for r in data_rows]
    elif args.NULL_missing_in:
        data_rows = [[None if v=='NULL' else v for v in r] for r in data_rows]
    # define the renamed input headers
    input_headers = Util.get_stripped_lines(raw_input_headers)
    if len(input_headers) < len(header):
        msg = 'each input header should be explicitly (re)named'
        raise ValueError(msg)
    if len(header) < len(input_headers):
        msg = 'more renamed headers than input headers'
        raise ValueError(msg)
    for h in input_headers:
        if not Carbone.is_valid_header(h):
            msg = 'invalid column header: %s' % h
            raise ValueError(msg)
    # force IC prefix for non-missing elements in the first column if requested
    if args.clean_isolates:
        data_rows = Carbone.clean_isolate_table(data_rows)
    # define the ordered output headers
    output_headers = Util.get_stripped_lines(raw_output_headers)
    bad_output_headers = set(output_headers) - set(input_headers)
    if bad_output_headers:
        msg_a = 'unrecognized output column headers: '
        msg_b = ', '.join(bad_output_headers)
        raise ValueError(msg_a + msg_b)
    # define the order of the output data columns
    h_to_i = dict((h, i) for i, h in enumerate(input_headers))
    # build the output data rows by reordering the columns
    data_rows = [[row[h_to_i[h]] for h in output_headers] for row in data_rows]
    # deal with missing data by skipping rows or replacing elements
    table = []
    for row in data_rows:
        if args.remove_missing_out and (None in row):
            continue
        elif args.NA_missing_out:
            row = ['NA' if x is None else x for x in row]
        table.append(row)
    # add row index labels for R compatibility if requested
    if args.add_indices:
        table = [[i+1] + row for i, row in enumerate(table)]
    # begin writing the R table
    out = StringIO()
    # write the table header
    print >> out, '\t'.join(output_headers)
    # write the table
    for row in table:
        print >> out, '\t'.join(str(x) for x in row)
    # return the table
    return out.getvalue()
