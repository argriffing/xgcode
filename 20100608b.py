"""Create an R table with fungus isolate info and with population PCA vectors.

The .hud file provides the names of the OTUs.
The amdS_PCA_Info.csv file provides other info.
Output is an R table suitable for plotting.
"""

from StringIO import StringIO
import sys
import os
import csv

import numpy as np
import argparse

from SnippetUtil import HandlingError
import Form
from Form import RadioItem
import Util
import EigUtil
import Carbone

g_default_hud_string = """
IC1 1 1 1 1
IC2 1 1 1 0
IC3 0 0 0 0
""".strip()

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

def do_pca(hud_lines):
    """
    @param hud_lines: lines of a .hud file
    @return: names, scaled vectors
    """
    # get the ordered names from the .hud file
    words = Carbone.get_words(hud_lines)
    names = [w.name for w in words]
    # create the floating point count matrix
    C_full = np.vstack([w.v for w in words])
    m_full, n_full = C_full.shape
    # remove invariant columns
    C = np.vstack([v for v in C_full.T if len(set(v))>1]).T
    # get the shape of the matrix
    m, n = C.shape
    # get the column means
    u = C.mean(axis=0)
    # get the centered and normalized counts matrix
    M = (C - u) / np.sqrt(u * (1 - u))
    # construct the sample covariance matrix
    X = np.dot(M, M.T) / n
    # get the eigendecomposition of the covariance matrix
    evals, evecs = EigUtil.eigh(X)
    # scale the eigenvectos by the eigenvalues
    pcs = [w*v for w, v in zip(evals, evecs)]
    return names, pcs

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
            Form.ContentDisposition()]
    return form_objects

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    text = process(fs.hud.splitlines(), fs.info.splitlines())
    disposition = "%s; filename=%s" % (fs.contentdisposition, 'data.table') 
    response_headers = [
            ('Content-Type', 'text/plain'),
            ('Content-Disposition', disposition)]
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

def process(hud_lines, info_lines):
    hud_lines = Util.get_stripped_lines(hud_lines)
    info_lines = Util.get_stripped_lines(info_lines)
    # read the .hud file and extract names and principal components
    names, pcs = do_pca(hud_lines)
    # check for sufficient number of eigenvectors
    if len(pcs) < 3:
        raise HandlingError('not enough principal components')
    # extract temperatures from the .csv file
    rows = list(csv.reader(info_lines))
    header, data_rows = rows[0], rows[1:]
    otu_to_info = {}
    for row in data_rows:
        otu = 'IC' + row[0]
        try:
            info = [
                    row_to_species(row),
                    row_to_location(row),
                    row_to_temperature(row),
                    row_to_precipitation(row)]
        except MissingError, e:
            continue
        otu_to_info[otu] = info
    # write the R table
    out = StringIO()
    #h = ('otu', 'species', 'location', 'temperature', 'precipitation',
        #'pc1', 'pc2', 'pc3', 'species.symbol', 'location.symbol')
    h = ('otu', 'species', 'location', 'temperature', 'precipitation',
        'pc1', 'pc2', 'pc3')
    print >> out, '\t'.join(h)
    for i, name in enumerate(names):
        if name in otu_to_info:
            info = otu_to_info[name]
            rowpcs = [pcs[0][i], pcs[1][i], pcs[2][i]]
            #species_symbol = ['AfL', 'Aa', 'Ac', 'Ano',
                    #'Ao', 'Ap', 'AfX', 'At',
                    #'Afu', 'Aso', 'AfS'].index(info[0]) + 1
            #symbols = [species_symbol]
            #symbols = [species_symbol, location_symbol]
            #location_symbol = ['AfL', 'Aa', 'Ac', 'Ano',
                    #'Ao', 'Ap', 'AfX', 'At',
                    #'Afu', 'Aso', 'AfS'].index(info[0]) + 1
            #row = [i+1, name] + info + rowpcs + symbols
            row = [i+1, name] + info + rowpcs
            print >> out, '\t'.join(str(x) for x in row)
    return out.getvalue()


def main(args):
    with open(os.path.expanduser(args.hud)) as fin_hud:
        with open(os.path.expanduser(args.csv)) as fin_info:
            print process(fin_hud, fin_info)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hud', help='.hud file')
    parser.add_argument('--info', help='a .csv like amdS_PCA_Info.csv')
    main(parser.parse_args())
