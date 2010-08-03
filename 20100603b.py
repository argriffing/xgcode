"""Create a fungus .ind file for location.

Create an .ind location file from a .hud and a amdS_PCA_Info.csv file.
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
import FormOut
import Util
import hud

g_tags = ['carbone_lab']

g_default_hud_string = """
IC1 1 1 1 0
IC2 1 1 1 0
IC3 1 0 1 0
""".strip()

g_default_info_lines = [
        '"IC","Haplo","Location","Temp (C)","Precip (mm)","Species",'
            '"B1","B2","G1","G2","OMST"',
        '"1","H42","NC","15","600","Ap","+","+","+","+","-"',
        '"2","H42","GA","30","700","Ap","+","+","+","+","-"',
        '"3","*","GA","45","800","Ap","+","+","+","+","-"']

g_default_info_string = '\n'.join(g_default_info_lines)

class DataRowError(Exception):
    def __init__(self, row, e):
        lines = ['Error in data row ' + str(row), str(e)]
        msg = '\n'.join(lines)
        Exception.__init__(self, msg)


def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('hud',
                'a list of OTUs names and binary character vectors',
                g_default_hud_string),
            Form.MultiLine('info',
                'amdS_PCA_Info.csv lines',
                g_default_info_string),
            Form.SingleLine('location',
                'control location',
                'GA'),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.EigenstratInd('location')

def get_response_content(fs):
    return process(
            fs.hud.splitlines(), fs.info.splitlines(), fs.location) + '\n'

def get_location_info(data_rows, control_location):
    """
    Asterisk is missing data.
    @param data_rows: rows of string elements
    @param control_location: control location string
    @return: location case and control OTU sets
    """
    cases = set()
    controls = set()
    for row in data_rows:
        try:
            otu = 'IC' + row[0]
            location = row[2]
            if location == '*':
                continue
            if location != control_location:
                cases.add(otu)
            else:
                controls.add(otu)
        except Exception, e:
            raise DataRowError(row, e)
    return cases, controls

def process(hud_lines, info_lines, location):
    """
    @param hud_lines: lines of a .hud file
    @param info_lines: lines of a phenotype .csv file
    @param location: the control location string
    """
    out = StringIO()
    # extract name order from the .hud file
    names, hud_data = hud.decode(hud_lines)
    # read the csv file
    rows = list(csv.reader(info_lines))
    header, data_rows = rows[0], rows[1:]
    cases, controls = get_location_info(data_rows, location)
    # write the .ind file contents
    for name in names:
        gender = 'U'
        classification = 'Ignore'
        if name in cases:
            classification = 'Case'
        elif name in controls:
            classification = 'Control'
        row = [name, gender, classification]
        print >> out, '\t'.join(row)
    return out.getvalue().rstrip()

def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        with open(os.path.expanduser(args.csv)) as fin_info:
            print process(fin_hud, fin_info, args.location)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', help='.hud file')
    parser.add_argument('--info', help='a .csv like amdS_PCA_Info.csv')
    parser.add_argument('--location', default='GA',
            help='control location')
    main(parser.parse_args())
