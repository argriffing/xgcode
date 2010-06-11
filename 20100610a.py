"""Plot some fungus info given an R table.
"""

from StringIO import StringIO
import os
import subprocess
from subprocess import PIPE
import tempfile

import argparse

from SnippetUtil import HandlingError
import Form
from Form import RadioItem
import Util
import iterutils

g_default_rows = [
        ['otu', 'species', 'location', 'temperature', 'precipitation',
            'pc1', 'pc2', 'pc3'],
        [1, 'IC100', 'Ap', 'GA', 15.0, 600.0,
            -2.8053476259, 0.556532380058, -6.17891756957],
        [2, 'IC101', 'Ap', 'GA', 15.0, 600.0,
            -2.8053476259, 0.556532380058, -6.17891756956],
        [3, 'IC102', 'Ap', 'GA', 15.0, 600.0,
            -2.80534762591, 0.556532380059, -6.17891756957],
        [455, 'IC577', 'Ac', 'AR', 25.0, 400.0,
            -13.7544736082, -7.16259232881, 7.0902951321],
        [456, 'IC580', 'Ac', 'AR', 25.0, 400.0,
            3.56768959361, 0.385873934264, 1.23981735331],
        [457, 'IC591', 'Ac', 'AR', 25.0, 400.0,
            -11.6455270418, -5.710582374, 5.60835091179]]

g_default_lines = ['\t'.join(str(x) for x in row) for row in g_default_rows]
g_default_string = '\n'.join(g_default_lines)

def my_mktemp():
    """
    The old tempfile.mktemp is deprecated.
    This function will create the name of a file
    that an external process can create.
    @return: the name of a nonexistent temporary file
    """
    f = tempfile.NamedTemporaryFile(delete=False)
    name = f.name
    f.close()
    os.unlink(name)
    return name

class MissingError(Exception): pass

def get_form():
    """
    @return: the body of a form
    """
    form_objects = [
            Form.MultiLine('table', 'R table', g_default_string),
            Form.Float('size', 'size of a plotted point', '1.5'),
            Form.SingleLine('legend_pos', 'position of the symbol legend',
                '0 0 0'),
            Form.RadioGroup('shape',
                'the variable that defines the shape of a plotted point', [
                    RadioItem('species', 'species', True),
                    RadioItem('location', 'location')]),
            Form.RadioGroup('color',
                'the variable that defines the color of a plotted point', [
                    RadioItem('temperature', 'temperature', True),
                    RadioItem('precipitation', 'precipitation')]),
            Form.RadioGroup('imageformat', 'output image format', [
                RadioItem('png', 'png', True),
                RadioItem('pdf', 'pdf'),
                RadioItem('postscript', 'postscript'),
                RadioItem('svg', 'svg')]),
            Form.ContentDisposition()]
    return form_objects

def create_filename(args):
    """
    @param args: user specified arguments from the web interface
    """
    fmt_to_ext = {
            'svg' : 'svg',
            'png' : 'png',
            'pdf' : 'pdf',
            'postscript' : 'ps'}
    ext = fmt_to_ext[args.imageformat]
    row = [args.shape, args.color, 'pca', '3d', ext]
    return '.'.join(row)

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    contents = process(fs, fs.table.splitlines())
    filename = create_filename(fs)
    format_to_content_type = {
            'svg':'image/svg+xml',
            'png':'image/png',
            'pdf':'application/pdf',
            'postscript':'application/postscript'}
    content_type = format_to_content_type[fs.imageformat]
    response_headers = [('Content-Type', content_type)]
    disposition = "%s; filename=%s" % (fs.contentdisposition, filename) 
    response_headers.append(('Content-Disposition', disposition)) 
    return response_headers, contents

def get_unique_species(table_lines):
    """
    @param table_lines: sequence of stripped lines of the R table
    @return: a list of unique species names
    """
    header, data_lines = table_lines[0], table_lines[1:]
    columns = zip(*[line.split() for line in data_lines])
    return list(iterutils.unique_everseen(columns[2]))

def get_unique_locations(table_lines):
    """
    @param table_lines: sequence of stripped lines of the R table
    @return: a list of unique location names
    """
    header, data_lines = table_lines[0], table_lines[1:]
    columns = zip(*[line.split() for line in data_lines])
    return list(iterutils.unique_everseen(columns[3]))

def get_augmented_table_lines(table_lines):
    header, data_lines = table_lines[0], table_lines[1:]
    n = len(data_lines)
    rows = [line.split() for line in data_lines]
    unique_species = get_unique_species(table_lines)
    unique_locations = get_unique_locations(table_lines)
    aug_header = '\t'.join([header, 'species.symbol', 'location.symbol'])
    aug_data_lines = []
    for row in rows:
        species = row[2]
        location = row[3]
        species_s = str(unique_species.index(species) + 1)
        location_s = str(unique_locations.index(location) + 1)
        aug_line = '\t'.join(row + [species_s, location_s])
        aug_data_lines.append(aug_line)
    return [aug_header] + aug_data_lines


def get_script(args, table_lines, temp_plot_filename, temp_table_filename):
    """
    @param args: from cmdline or web
    @param table_lines: sequence of stripped lines of the R table
    @param temp_plot_name: a pathname
    @param temp_table_name: a pathname
    """
    # get the symbol legend location
    legend_pos = tuple(float(x) for x in args.legend_pos.split())
    # get the unique locations and species
    unique_species = get_unique_species(table_lines)
    unique_locations = get_unique_locations(table_lines)
    species_string = ', '.join("'%s'" % x for x in unique_species)
    locations_string = ', '.join("'%s'" % x for x in unique_locations)
    if args.shape == 'species':
        symbol_legend_string = species_string
    elif args.shape == 'location':
        symbol_legend_string = locations_string
    else:
        raise ValueError('invalid args.shape')
    if args.color == 'temperature':
        color_legend_string = 'celsius temperature'
    elif args.color == 'precipitation':
        color_legend_string = 'millimeters of precipitation'
    else:
        raise ValueError('invalid args.color')
    rcodes = [
        "require('scatterplot3d')",
        "mytable <- read.table('%s')" % temp_table_filename,
        "%s('%s')" % (args.imageformat, temp_plot_filename),
        # rename some variables for compatibility with the template
        "Year <- mytable$pc1",
        "Latitude <- mytable$pc2",
        "Risk <- mytable$pc3",
        "Prec <- mytable$%s" % args.color,
        # stack two plots vertically
        "layout(cbind(1:2, 1:2), heights = c(7, 1))",
        # create the color gradient
        "prc <- hsv((prc <-0.7 * Prec / diff(range(Prec))) - min(prc) + 0.3)",
        # create the scatterplot
        "s3d <- scatterplot3d(Year, Latitude, Risk, mar = c(5, 3, 4, 3),",
        "xlab = 'PC1', ylab = 'PC2', zlab = 'PC3',",
        "type = 'p', pch = ' ')",
        # define symbols colors and sizes
        "s3d$points(Year, Latitude, Risk,",
        "pch=mytable$%s.symbol," % args.shape,
        "bg=prc, col=prc, cex=%s)" % args.size,
        # define x y and z as Year, Latitude and Risk
        "s3d.coords <- s3d$xyz.convert(Year, Latitude, Risk)",
        # symbol legend
        "legend(s3d$xyz.convert(%s, %s, %s)," % legend_pos,
        "pch=1:%s, yjust = 0," % len(symbol_legend_string),
        "legend=c(%s)," % symbol_legend_string,
        "cex = 1.1)",
        # set margins
        "par(mar=c(5, 3, 0, 3))",
        # draw the plot
        "plot(seq(min(Prec), max(Prec), length = 100),",
        "rep(0, 100), pch = 15,",
        "axes = FALSE,",
        "xlab = '%s'," % color_legend_string,
        "ylab = '', col = hsv(seq(0.3, 1, length = 100)))",
        # draw the color legend
        "axis(1)",
        # write the plot
        "dev.off()"]
    return '\n'.join(rcodes)

def process(args, table_lines):
    """
    @param args: command line or web input
    @param table_lines: input lines
    @return: the image data as a string
    """
    table_lines = Util.get_stripped_lines(table_lines)
    # Create a temporary data table file for R.
    f_temp_table = tempfile.NamedTemporaryFile(delete=False)
    augmented_table_lines = get_augmented_table_lines(table_lines)
    print >> f_temp_table, '\n'.join(augmented_table_lines)
    f_temp_table.close()
    # Create a temporary pathname for the plot created by R.
    temp_plot_name = my_mktemp()
    # Create a temporary R script file.
    f_temp_script = tempfile.NamedTemporaryFile(delete=False)
    script = get_script(args, table_lines, temp_plot_name, f_temp_table.name)
    print >> f_temp_script, script
    f_temp_script.close()
    # Call R.
    cmd = ['R', 'CMD', 'BATCH', '--vanilla', f_temp_script.name]
    proc = subprocess.Popen(cmd, stdout=PIPE, stderr=PIPE)
    r_out, r_err = proc.communicate()
    # Delete the temporary data table file.
    os.unlink(f_temp_table.name)
    # Delete the temporary script file.
    os.unlink(f_temp_script.name)
    # Read the image file.
    try:
        with open(temp_plot_name, 'rb') as fin:
            image_data = fin.read()
    except IOError, e:
        raise HandlingError('the R call seems to not have created the plot')
    # Delete the temporary image file.
    os.unlink(temp_plot_name)
    # Return the image data as a string.
    return image_data

def main(args):
    with open(os.path.expanduser(args.table)) as fin:
        sys.stdout.write(process(args, fin))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', required=True,
            help='an R table with a specific format')
    parser.add_argument('--imageformat', default='png',
            choices=['png', 'pdf', 'postscript', 'svg'],
            help='output image file format')
    parser.add_argument('--shape', required=True,
            choices=['species', 'location'],
            help='the variable that defines the shape of a plotted point')
    parser.add_argument('--color', required=True,
            choices=['temperature', 'precipitation'],
            help='the variable that defines the color of a plotted point')
    parser.add_argument('--size', default=1.5, type=float,
            help='the size of each plotted point')
    parser.add_argument('--legend_pos', default='0 0 0',
            help='position of the symbol legend')
    main(parser.parse_args())
