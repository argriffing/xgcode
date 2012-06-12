"""
Modify a sequence alignment within a BEAST xml, by taking a column subset.

The number of MCMC steps is also set.
The web interface modifies the primates.xml from the BEAST tutorial.
"""

from StringIO import StringIO
import argparse
import sys

from lxml import etree

import Form
import FormOut
import const

g_primates_xml = const.read('20120605a')


def get_form():
    form_objects = [
            Form.IntegerInterval(
                'start_pos', 'stop_pos', 'alignment interval',
                1, 4263, low=1, high=4263),
            Form.Integer('nsteps', 'MCMC steps', 8000, low=1),
            Form.SingleLine('alignment_id', 'alignment id', 'alignment'),
            Form.SingleLine('mcmc_id', 'mcmc id', 'mcmc'),
            ]
    return form_objects

def get_form_out():
    return FormOut.BeastXml()

def modify_taxon_sequence(taxon_element, start_pos, stop_pos):
    sequence = taxon_element.tail.strip()
    taxon_element.tail = sequence[start_pos-1 : stop_pos]

def process(fs, xmldata):
    """
    @param fs: a fieldstorage-like object
    @param xmldata: contents of an xml file as a multiline string
    @return: contents of the new xml file as a multiline string
    """
    # get user input
    start_pos = fs.start_pos
    stop_pos = fs.stop_pos
    nsteps = fs.nsteps
    alignment_id = fs.alignment_id
    # read the xml tree
    tree = etree.parse(StringIO(xmldata))
    # modify the number of mcmc steps
    for event, element in etree.iterwalk(tree, tag='mcmc'):
        if element.get('id') == fs.mcmc_id:
            element.set('chainLength', str(fs.nsteps))
    # modify the sequences within the alignment
    for event, element in etree.iterwalk(tree, tag='alignment'):
        if element.get('id') == fs.alignment_id:
            for seq_element in element:
                if seq_element.tag != 'sequence':
                    continue
                for taxon_element in seq_element:
                    if taxon_element.tag != 'taxon':
                        continue
                    modify_taxon_sequence(
                            taxon_element, start_pos, stop_pos)
    # write the xml tree
    return etree.tostring(tree)

def get_response_content(fs):
    return process(fs, g_primates_xml)

def main(args):
    if args.infile is None:
        xmldata_in = sys.stdin.read()
    else:
        with open(args.infile) as fin:
            xmldata_in = fin.read()
    xmldata_out = process(args, xmldata_in)
    if args.outfile is None:
        sys.stdout.write(xmldata_out)
    else:
        with open(args.outfile) as fout:
            fout.write(xmldata_out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', help='input xml file')
    parser.add_argument('-o', '--outfile', help='output xml file')
    parser.add_argument('--start_pos',
            default=1, type=int,
            help='alignment interval start position')
    parser.add_argument('--stop_pos',
            default=1, type=int,
            help='alignment interval stop position')
    parser.add_argument('--nsteps',
            default=8000, type=int,
            help='number of MCMC steps')
    parser.add_argument('--alignment_id',
            default='alignment',
            help='alignment id')
    parser.add_argument('--mcmc_id',
            default='mcmc',
            help='mcmc id')
    main(parser.parse_args())

