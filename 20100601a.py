"""Create a fungus precipitation .ind file
Create an .ind precipitation file from a .hud and a amdS_PCA_Info.csv file.

The .hud file provides the names of the OTUs.
The amdS_PCA_Info.csv file provides the 'case-control' status,
representing binarized location, temperature, or precipitation.
The output file is in eigenstrat format.
"""

from StringIO import StringIO
import sys
import os
import csv

import argparse

from SnippetUtil import HandlingError
import Form
import Util
import Carbone

g_default_hud_string = """
foo 1 1 1
bar 1 1 1
baz 1 0 1
""".strip()

g_default_info_lines = [
        '"IC","Haplo","Location","Temp (C)","Precip (mm)","Species",'
            '"B1","B2","G1","G2","OMST"',
        '"1","H42","GA","15","600","Ap","+","+","+","+","-"',
        '"2","H42","GA","15","600","Ap","+","+","+","+","-"',
        '"3","*","GA","15","600","Ap","+","+","+","+","-"']

g_default_info_string = '\n'.join(g_default_info_lines)

def get_validated_words(lines):
    lines = Util.get_stripped_lines(lines)
    words = [Carbone.Word(line) for line in lines]
    Carbone.validate_words(words)
    return words


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'a list of OTUs names and binary character vectors',
                g_default_hud_string),
            Form.Multiline('info',
                'amdS_PCA_Info.csv lines',
                g_default_info_lines),
            Form.Float('threshold',
                    'precipitation threshold (mm)',
                    '750.0')]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    lines = Util.get_stripped_lines(StringIO(fs.hud))
    words = get_validated_words(lines)
    return [('Content-Type', 'text/plain')], text

def main(args):
    # extract names from the .hud file
    with open(args.hud) as fin:
        lines = Util.get_stripped_lines(fin)
    words = get_validated_words(lines)
    names = [word.name for word in words]
    with open(args.info) as fin:
        pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', help='.hud file')
    parser.add_argument('--info', help='a .csv like amdS_PCA_Info.csv')
    parser.add_argument('--threshold', type=float, default=750,
            help='precipitation threshold (mm)')
    args = parser.parse_args()
    main(args)




class DataRowError(Exception):
    def __init__(self, row, e):
        lines = ['Error in data row ' + str(row), str(e)]
        msg = '\n'.join(lines)
        Exception.__init__(self, msg)

def main(args):
    # get the names from the .hud file
    names = []
    with open(os.path.expanduser(args.hud)) as fin_hud:
        for line in util.gen_nonempty_stripped_lines(fin_hud):
            name, rest = line.split(None, 1)
            names.append(name)
    # open the csv file
    with open(os.path.expanduser(args.csv)) as fin_csv:
        # start reading the csv file
        rows = list(csv.reader(fin_csv))
        header, data_rows = rows[0], rows[1:]
        # get case and control OTU sets
        if args.environment == 'precipitation':
            cases, controls = get_precipitation_info(data_rows,
                    args.precipitation_threshold)
        elif args.environment == 'temperature':
            cases, controls = get_temperature_info(data_rows,
                    args.temperature_threshold)
        elif args.environment == 'location':
            cases, controls = get_location_info(data_rows,
                    args.control_location)
        else:
            msg = 'unrecognized environmental variable: ' + args.environment
            raise Exception(msg)
    # write the .ind file contents
    for name in names:
        gender = 'U'
        classification = 'Ignore'
        if name in cases:
            classification = 'Case'
        elif name in controls:
            classification = 'Control'
        row = [name, gender, classification]
        print '\t'.join(row)

def get_precipitation_info(data_rows, threshold):
    """
    Asterisk is missing data.
    @param data_rows: rows of string elements
    @param threshold: precipitation threshold in millimeters
    @return: precipitation case and control OTU sets
    """
    cases = set()
    controls = set()
    for row in data_rows:
        try:
            otu = 'IC' + row[0]
            precipitation = row[4]
            if precipitation == '*':
                continue
            if float(precipitation) < threshold:
                cases.add(otu)
            else:
                controls.add(otu)
        except Exception, e:
            raise DataRowError(row, e)
    return cases, controls


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precipitation_threshold', type=float,
            default=750,
            help='threshold classifying the amount of precipitation (mm)')
    parser.add_argument('--hud', required=True,
            help='an input .hud file')
    parser.add_argument('--csv', required=True,
            help = 'an input amdS_PCA_Info.csv file')
    args = parser.parse_args()
    main(args)
