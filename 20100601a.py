"""
Convert a .hud file and a amdS_PCA_Info.csv file to an .ind file.
The .hud file provides the names of the OTUs.
The amdS_PCA_Info.csv file provides the 'case-control' status,
representing binarized location, temperature, or precipitation.
The output file is in eigenstrat format.
"""

import os
import sys
import csv

import argparse

import util

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

def get_temperature_info(data_rows, threshold):
    """
    Asterisk is missing data.
    @param data_rows: rows of string elements
    @param threshold: temperature threshold in Celcius
    @return: temperature case and control OTU sets
    """
    cases = set()
    controls = set()
    for row in data_rows:
        try:
            otu = 'IC' + row[0]
            temperature = row[3]
            if temperature == '*':
                continue
            if float(temperature) < threshold:
                cases.add(otu)
            else:
                controls.add(otu)
        except Exception, e:
            raise DataRowError(row, e)
    return cases, controls

def get_location_info(data_rows, control_location):
    """
    Asterisk is missing data.
    @param data_rows: rows of string elements
    @param control_location: the location treated as the control
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
            if location == control_location:
                controls.add(otu)
            else:
                cases.add(otu)
        except Exception, e:
            raise DataRowError(row, e)
    return cases, controls


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--precipitation_threshold', type=float,
            default=750,
            help='threshold classifying the amount of precipitation (mm)')
    parser.add_argument('--temperature_threshold', type=float,
            default=22,
            help='threshold to classify the temperature (C)')
    parser.add_argument('--control_location',
            default='GA',
            help='use this location as the control location')
    parser.add_argument('--environment',
            choices=['location', 'temperature', 'precipitation'],
            default='temperature',
            help='the environmental variable of interest')
    parser.add_argument('--hud', required=True,
            help='an input .hud file')
    parser.add_argument('--csv', required=True,
            help = 'an input amdS_PCA_Info.csv file')
    args = parser.parse_args()
    main(args)
