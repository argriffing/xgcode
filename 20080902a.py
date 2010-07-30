"""Visualize MAPP p-values for one missense mutation per alignment column.

MAPP stands for Multivariate Analysis of Protein Polymorphism.
"""

import math
from StringIO import StringIO

import cairo

from SnippetUtil import HandlingError
import SnippetUtil
import Newick
import CairoUtil
import Codon
import Util
import Form
import FormOut
import const

g_mapp_output = const.read('20100730v')


class ColumnDataError(Exception):
    pass


class ColumnData:
    """
    Define data associated with the column of an alignment.
    """

    def __init__(self, gene, offset, wild, mutant):
        """
        @param gene: the name of a gene
        @param offset: the integer offset indexed starting at 1
        @param wild: the wild type amino acid
        @param mutant: the mutant type amino acid
        """
        if type(offset) is not type(1):
            raise ColumnDataError('invalid offset')
        if offset < 1:
            raise ColumnDataError('the offset must be greater or equal to one')
        if wild not in Codon.g_aa_letters:
            raise ColumnDataError('invalid wild type')
        if mutant not in Codon.g_aa_letters:
            raise ColumnDataError('invalid mutant type')
        self.gene = gene
        self.offset = offset
        self.wild = wild
        self.mutant = mutant

    def __str__(self):
        return self.gene + ':' + self.wild + str(self.offset) + self.mutant


def get_image_string(
        pvalue_lists, header_list, wild_list, mutant_list, image_format):
    """
    @param pvalue_lists: for each amino acid column, a list of aa pvalues
    @param wild_list: for each amino acid column, the wild type amino acid
    @param mutant_list: for each amino acid column, the mutant amino acid
    @param image_format: something like 'png'
    """
    # start with unrealistic image dimensions to get some font size information
    initial_width = 100
    initial_height = 100
    initial_cairo_helper = CairoUtil.CairoHelper(image_format)
    initial_surface = initial_cairo_helper.create_surface(
            initial_width, initial_height)
    initial_context = cairo.Context(initial_surface)
    initial_context.select_font_face(
            'monospace', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    initial_context.set_font_size(12)
    extents = initial_context.font_extents()
    ascent, descent, font_height, max_x_advance, max_y_advance = extents
    # get the blank image because multiple surfaces cannot exist simultaneously
    dummy_string = initial_cairo_helper.get_image_string()
    # use a standard image width
    image_width = 640
    # Use one row of text for each alignment column header,
    # and three rows for the axis.
    image_height = font_height * (3 + len(header_list))
    # create the context with realistic image dimensions
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(image_width, image_height)
    context = cairo.Context(surface)
    # set the font size to be used throughout the visualization
    context.select_font_face(
            'monospace', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    context.set_font_size(12)
    # draw a light gray background so you can see it in a white canvas
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # get the x advance of the longest header
    max_header_advance = max(
            context.text_extents(header)[4] for header in header_list)
    # draw the right-justified headers
    for i, header in enumerate(header_list):
        extents = context.text_extents(header)
        x_bearing, y_bearing, text_w, text_h, x_advance, y_advance = extents
        x_offset = max_header_advance - x_advance
        y_offset = font_height * (i + 1)
        context.move_to(x_offset, y_offset)
        context.show_text(header)
        context.stroke()
    # get the p-value range
    pvalue_min = min(min(pvalue_list) for pvalue_list in pvalue_lists)
    pvalue_max = max(max(pvalue_list) for pvalue_list in pvalue_lists)
    # get the pixel offset to start the x axis
    x_axis_offset = max_header_advance + max_x_advance
    # get the pixel width of the x axis
    x_axis_width = (image_width - max_x_advance) - x_axis_offset
    # get the number of pixels per log p-value
    pixels_per_log_pvalue = (
            x_axis_width / (math.log(pvalue_min) / math.log(10)))
    # draw the amino acid letters in gray,
    # but do not draw the wild type or the mutant amino acids
    context.save()
    context.set_source_rgb(.7, .7, .7)
    for header_index, header in enumerate(header_list):
        y_offset = font_height * (header_index + 1)
        for aa_index, aa_letter in enumerate(Codon.g_aa_letters):
            if aa_letter == wild_list[header_index]:
                continue
            if aa_letter == mutant_list[header_index]:
                continue
            pvalue = pvalue_lists[header_index][aa_index]
            log_pvalue = math.log(pvalue) / math.log(10)
            x_offset = x_axis_offset + pixels_per_log_pvalue * log_pvalue
            context.move_to(x_offset, y_offset)
            context.show_text(aa_letter)
            context.stroke()
    context.restore()
    # draw the x axis
    x_axis_y_offset = font_height * (len(header_list) + 0.5)
    context.save()
    context.set_line_cap(cairo.LINE_CAP_ROUND)
    context.move_to(x_axis_offset, x_axis_y_offset)
    context.line_to(x_axis_offset + x_axis_width, x_axis_y_offset)
    context.stroke()
    context.restore()
    # draw the x axis ticks
    context.save()
    context.set_line_cap(cairo.LINE_CAP_ROUND)
    log_pvalue = 0
    while math.exp(log_pvalue) > pvalue_min:
        x_offset = x_axis_offset + log_pvalue * pixels_per_log_pvalue
        context.move_to(x_offset, x_axis_y_offset)
        context.line_to(x_offset, x_axis_y_offset + font_height * 0.5)
        context.stroke()
        log_pvalue -= 1
    context.restore()
    # draw the x axis tick labels
    log_pvalue = 0
    while math.exp(log_pvalue) > pvalue_min:
        x_offset = x_axis_offset + log_pvalue * pixels_per_log_pvalue
        y_offset = font_height * (len(header_list) + 2)
        context.move_to(x_offset, y_offset)
        context.show_text('1E%d' % log_pvalue)
        context.stroke()
        log_pvalue -= 1
    # make a utility dictionary
    aa_letter_to_index = {}
    for i, aa_letter in enumerate(Codon.g_sorted_aa_letters):
        aa_letter_to_index[aa_letter] = i
    # draw the wild type amino acids in blue
    context.save()
    context.set_source_rgb(0.2, 0.2, 1.0)
    for header_index, aa_letter in enumerate(wild_list):
        y_offset = font_height * (header_index + 1)
        pvalue = pvalue_lists[header_index][aa_letter_to_index[aa_letter]]
        log_pvalue = math.log(pvalue) / math.log(10)
        x_offset = x_axis_offset + pixels_per_log_pvalue * log_pvalue
        context.move_to(x_offset, y_offset)
        context.show_text(aa_letter)
        context.stroke()
    context.restore()
    # draw the mutant amino acids in red
    context.save()
    context.set_source_rgb(1.0, 0.2, 0.2)
    for header_index, aa_letter in enumerate(mutant_list):
        y_offset = font_height * (header_index + 1)
        pvalue = pvalue_lists[header_index][aa_letter_to_index[aa_letter]]
        log_pvalue = math.log(pvalue) / math.log(10)
        x_offset = x_axis_offset + pixels_per_log_pvalue * log_pvalue
        context.move_to(x_offset, y_offset)
        context.show_text(aa_letter)
        context.stroke()
    context.restore()
    # get the image string
    return cairo_helper.get_image_string()

def get_form():
    """
    @return: the body of a form
    """
    # define some default data
    column_info_list = (
            ColumnData('TMEM195', 279, 'F', 'L'),
            ColumnData('SHANK3', 552, 'R', 'W'),
            ColumnData('ARHGAP6', 76, 'G', 'D'),
            ColumnData('GRPR', 261, 'R', 'L'),
            ColumnData('DRP2', 203, 'T', 'M'),
            ColumnData('PLXNB3', 981, 'R', 'H'))
    data_lines = [str(info) for info in column_info_list]
    # define the default wild and mutant amino acids
    wild_string = ''.join(info.wild for info in column_info_list)
    mutant_string = ''.join(info.mutant for info in column_info_list)
    # define the form objects
    form_objects = [
            Form.MultiLine('mapp', 'MAPP output',
                g_mapp_output.strip()),
            Form.MultiLine('headers', 'alignment column headers',
                '\n'.join(data_lines)),
            Form.SingleLine('wild', 'wild type amino acids', wild_string),
            Form.SingleLine('mutant', 'mutant amino acids', mutant_string),
            Form.ImageFormat(),
            Form.ContentDisposition()]
    return form_objects

def get_form_out():
    return FormOut.Image('out', [])

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response type
    response_headers = []
    # read the headers
    headers = Util.get_stripped_lines(fs.headers.splitlines())
    # read the wild type amino acid list
    wild_list = list(fs.wild.strip())
    extra_letters = set(wild_list) - set(Codon.g_aa_letters)
    if extra_letters:
        msg_a = 'some invalid wild type amino acids were specified: '
        msg_b = str(tuple(extra_letters))
        raise HandlingError(msg_a + msg_b)
    # read the mutant amino acid list
    mutant_list = list(fs.mutant.strip())
    extra_letters = set(mutant_list) - set(Codon.g_aa_letters)
    if extra_letters:
        msg_a = 'some invalid mutant amino acids were specified: '
        msg_b = str(tuple(extra_letters))
        raise HandlingError(msg_a + msg_b)
    # read the tab separated MAPP output
    tsv_lists = []
    for line in StringIO(fs.mapp):
        if line.strip():
            tsv_list = [element.strip() for element in line.split('\t')]
            tsv_lists.append(tsv_list)
    # check input consistency
    if len(headers) != len(tsv_lists) - 1:
        msg_a = 'the number of headers should be '
        msg_b = 'one fewer than the number of MAPP lines'
        raise HandlingError(msg_a + msg_b)
    if len(mutant_list) != len(tsv_lists) - 1:
        msg_a = 'the number of mutant amino acids should be '
        msg_b = 'one fewer than the number of MAPP lines'
        raise HandlingError(msg_a + msg_b)
    if len(wild_list) != len(tsv_lists) - 1:
        msg_a = 'the number of wild type amino acids should be '
        msg_b = 'one fewer than the number of MAPP lines'
        raise HandlingError(msg_a + msg_b)
    length_set = set(len(tsv_list) for tsv_list in tsv_lists)
    if length_set != set([54]):
        msg_a = 'each line in the MAPP output should have '
        msg_b = '54 tab separated values'
        raise HandlingError(msg_a + msg_b)
    # read the p-values
    pvalue_lists = []
    for tsv_list in tsv_lists[1:]:
        pvalue_list = []
        for element in tsv_list[32:52]:
            invalid_pvalue_error = HandlingError('invalid p-value: ' + element)
            try:
                pvalue = float(element)
            except ValueError:
                raise invalid_pvalue_error
            if pvalue <= 0.0 or pvalue > 1.0:
                raise invalid_pvalue_error
            pvalue_list.append(pvalue)
        pvalue_lists.append(pvalue_list)
    # get some options
    ext = Form.g_imageformat_to_ext[fs.imageformat]
    filename = 'mapp.' + ext
    contenttype = Form.g_imageformat_to_contenttype[fs.imageformat]
    contentdisposition = '%s; filename=%s' % (fs.contentdisposition, filename)
    # draw the image
    try:
        image_string = get_image_string(
                pvalue_lists, headers, wild_list, mutant_list, ext)
    except CairoUtil.CairoUtilError, e:
        raise HandlingError(e)
    # return the response
    response_headers = [
            ('Content-Type', contenttype),
            ('Content-Disposition', contentdisposition)]
    return response_headers, image_string
