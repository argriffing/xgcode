"""Fungus csv janitor.

Clean the fungus csv file using ad hoc codes.
The output is an R table.
"""

from StringIO import StringIO
import csv
import re

from SnippetUtil import HandlingError
import Form
import Util
import Carbone

g_info_lines = [
        '"IC","Haplo","Location","Temp (C)","Precip (mm)","Species",'
            '"B1","B2","G1","G2","OMST"',
        '"1","H42","GA","15","600","Ap","+","+","+","+","-"',
        '"2","H42","GA","30","700","Ap","+","+","+","+","-"',
        '"3","*","GA","45","800","Ap","+","+","+","+","-"']

g_input_headers = [
        'otu', 'haplotype', 'location', 'temperature',
        'precipitation', 'species',
        'b1', 'b2', 'g1', 'g2', 'omst']

g_output_headers = [
        'temperature', 'precipitation', 'otu']

g_column_header_pattern = r'^[a-zA-Z][a-zA-Z0-9]*$'


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('info', 'amdS_PCA_Info.csv lines',
                '\n'.join(g_info_lines)),
            Form.MultiLine('input_headers', 'renamed column headers',
                '\n'.join(g_input_headers)),
            Form.MultiLine('output_headers', 'ordered output headers',
                '\n'.join(g_output_headers)),
            Form.RadioGroup('mdefgroup', 'missing input data', [
                Form.RadioItem('star_missing_in',
                    '* is interpreted as missing data', True),
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

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs, fs.info.splitlines(),
            fs.input_headers.splitlines(), fs.output_headers.splitlines())
    filename = 'fungus.table'
    response_headers = [('Content-Type', 'text/plain')] 
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers.append(('Content-Disposition', disposition)) 
    return response_headers, text


def process(args, raw_info_lines, raw_input_headers, raw_output_headers):
    info_lines = Util.get_stripped_lines(raw_info_lines)
    # extract info from the .csv file
    rows = list(csv.reader(info_lines))
    # the number of columns should be consistent among rows
    if len(set(len(row) for row in rows)) != 1:
        msg = 'the number of columns should be consistent among rows'
        raise ValueError(msg)
    # break the list of rows into a header row and data rows
    header, data_rows = rows[0], rows[1:]
    # account for missing input data
    if args.star_missing_in:
        data_rows = [[None if v=='*' else v for v in r] for r in data_rows]
    # define the renamed input headers
    input_headers = Util.get_stripped_lines(raw_input_headers)
    if len(input_headers) < len(header):
        msg = 'each input header should be explicitly (re)named'
        raise ValueError(msg)
    if len(header) < len(input_headers):
        msg = 'more renamed headers than input headers'
        raise ValueError(msg)
    for h in input_headers:
        if not re.match(g_column_header_pattern, h):
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
