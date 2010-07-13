"""Show unique fungus locations and species.

The amdS_PCA_Info.csv file provides miscellaneous fungus info.
"""

from StringIO import StringIO
import os
import csv

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
from Form import RadioItem
import Util

g_tags = ['carbone_lab']

g_default_info_lines = [
        '"IC","Haplo","Location","Temp (C)","Precip (mm)","Species",'
            '"B1","B2","G1","G2","OMST"',
        '"1","H42","GA","15","600","Ap","+","+","+","+","-"',
        '"2","H42","GA","30","700","Ap","+","+","+","+","-"',
        '"3","*","GA","45","800","Ap","+","+","+","+","-"']

g_default_info_string = '\n'.join(g_default_info_lines)

class MissingError(Exception): pass

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
            Form.MultiLine('info',
                'amdS_PCA_Info.csv lines',
                g_default_info_string)]
    return form_objects

def get_form_out():
    return FormOut.Report('out.txt', [])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    response_headers = [('Content-Type', 'text/plain')] 
    text = process(fs.info.splitlines())
    return response_headers, text

def row_to_temperature(row):
    t = row[3]
    if t == '*':
        raise MissingError
    try:
        temperature = float(t)
    except ValueError, e:
        raise DataRowError(row, e)
    return temperature

def row_to_precipitation(row):
    p = row[4]
    if p == '*':
        raise MissingError
    try:
        precipitation = float(p)
    except ValueError, e:
        raise DataRowError(row, e)
    return precipitation

def row_to_species(row):
    species = row[5]
    if species == '*':
        raise MissingError
    return species

def row_to_location(row):
    location = row[2]
    if location == '*':
        raise MissingError
    return location

def process(info_lines):
    info_lines = Util.get_stripped_lines(info_lines)
    # init the species and location set
    species_set = set()
    location_set = set()
    # extract info from the .csv file
    rows = list(csv.reader(info_lines))
    header, data_rows = rows[0], rows[1:]
    for row in data_rows:
        try:
            info = [
                    row_to_species(row),
                    row_to_location(row),
                    row_to_temperature(row),
                    row_to_precipitation(row)]
        except MissingError, e:
            continue
        species_set.add(row_to_species(row))
        location_set.add(row_to_location(row))
    # write a summary
    out = StringIO()
    print >> out, 'found', len(species_set), 'unique species:'
    for s in species_set:
        print >> out, s
    print >> out
    print >> out, 'found', len(location_set), 'unique locations:'
    for s in location_set:
        print >> out, s
    return out.getvalue()


def main(args):
    with open(os.path.expanduser(args.info)) as fin_info:
        print process(fin_hud, fin_info)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--info', help='a .csv like amdS_PCA_Info.csv')
    main(parser.parse_args())
