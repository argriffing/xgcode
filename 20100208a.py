"""Split a file by the first whitespace separated element on each line. [NOWEB]

This is not so useful as a web app.
The file is split by lines into multiple files.
"""


from StringIO import StringIO
import os

import argparse

from SnippetUtil import HandlingError
import Form
import FormOut
import iterutils


def get_form():
    """
    @return: the body of a form
    """
    return []

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    return 'this is not web accessible\n'

def main(args):
    # read the arguments
    input_filename = os.path.abspath(os.path.expanduser(args.infile))
    output_directory = os.path.abspath(os.path.expanduser(args.outdir))
    force = args.force
    # make sure that the output directory exists
    if not os.path.isdir(output_directory):
        if force:
            os.makedirs(output_directory)
    if not os.path.isdir(output_directory):
        msg = 'output directory does not exist: ' + output_directory
        raise Exception(msg)
    # Open files and write as the filenames are discovered
    # unless this is a dry run.
    name_to_path = {}
    name_to_fout = {}
    with open(input_filename) as fin:
        for line in fin:
            short_row = line.split(None, 1)
            if not short_row:
                continue
            # get the name
            name = short_row[0]
            # do new name stuff if the name looks new
            if name not in name_to_path:
                # create the pathname
                output_filename = args.out_prefix + name + args.out_suffix
                fpath = os.path.join(output_directory, output_filename)
                name_to_path[name] = fpath
                # check for existence
                if not args.force:
                    if os.path.exists(fpath):
                        raise Exception('output file already exists: ' + fpath)
                # open for writing unless it is a dry run
                if not args.dryrun:
                    fout = open(fpath, 'w')
                    name_to_fout[name] = fout
            # If the open file exists then write the line.
            if name in name_to_fout:
                name_to_fout[name].write(line)
    # close any open files
    for fout in name_to_fout.values():
        fout.close()
    # if it is a dry run then list the pathnames
    if args.dryrun:
        for fpath in sorted(name_to_path.values()):
            print fpath


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('--force', action='store_true',
            help='overwrite existing files')
    parser.add_argument('--dryrun', action='store_true',
            help='list the names of files to be created')
    parser.add_argument('--outdir', default=os.getcwd(),
            help='write the files to this directory')
    parser.add_argument('--out_prefix', default='chromosome.',
            help='prefix added to the value name in the output filename')
    parser.add_argument('--out_suffix', default='.txt',
            help='suffix added to the value name in the output filename')
    args = parser.parse_args()
    main(args)
