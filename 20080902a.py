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

# sample MAPP output
g_mapp_output = """
Position	Column Score	Column p-value	Alignment	Gap Weight	Over Gap Weight Threshold	Hydropathy	Polarity	Charge	Volume	Free energy alpha	Free energy beta	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	A	C	D	E	F	G	H	I	K	L	M	N	P	Q	R	S	T	V	W	Y	Good Amino Acids	Bad Amino Acids
1	4.16	3.314E-1	'FFLVIFHVVAV-	0.1489702059588082	N	9.802E-1	9.481E-1	10E-1	10E-1	9.989E-1	9.927E-1	4.98	2.73	11.81	9.34	2.68	4.87	3.89	1.31	8.55	3.09	4.41	4.83	12.41	3.13	10.82	3.18	3.50	1.78	4.55	3.77	1.968E-1	7.09E-1	4.327E-3	1.439E-2	7.222E-1	2.114E-1	3.888E-1	9.836E-1	2.206E-2	6.03E-1	2.84E-1	2.165E-1	3.317E-3	5.916E-1	6.838E-3	5.77E-1	4.881E-1	9.346E-1	2.596E-1	4.172E-1	ACEFGHIKLMNQSTVWY	DPR
2	13.64	1.989E-3	'RRRR--RLR-RR	0.2370525894821035	N	9.963E-1	9.52E-1	8.222E-1	1.383E-1	1.189E-1	10E-1	18.51	21.62	13.66	14.32	10.53	28.25	10.14	10.23	2.67	4.35	5.23	13.41	39.30	5.94	0.78	21.21	19.73	17.57	14.88	13.62	3.604E-4	1.48E-4	1.97E-3	1.521E-3	7.879E-3	3.116E-5	9.525E-3	9.106E-3	7.245E-1	2.947E-1	1.674E-1	2.186E-3	4.434E-6	1.063E-1	9.989E-1	1.653E-4	2.503E-4	4.845E-4	1.231E-3	2.009E-3	KLMQR	ACDEFGHINPSTVWY
3	3.94	3.781E-1	'GGGGGG----TE	0.4380123975204958	N	9.998E-1	10E-1	9.996E-1	9.833E-1	10E-1	9.994E-1	5.05	4.73	3.21	3.28	5.22	1.95	3.32	7.28	3.64	6.71	3.76	2.47	5.69	2.86	4.18	0.89	2.45	7.18	4.11	4.12	1.881E-1	2.309E-1	5.673E-1	5.468E-1	1.689E-1	9.046E-1	5.364E-1	4.609E-2	4.504E-1	6.53E-2	4.205E-1	7.826E-1	1.247E-1	6.691E-1	3.279E-1	9.979E-1	7.866E-1	4.904E-2	3.406E-1	3.392E-1	ACDEFGHIKLMNPQRSTVWY	
4	22.61	1.143E-4	'RRRRR-RRRRRR	0.06568686262747449	N	2.95E-2	1.623E-2	1.734E-2	9.106E-2	1.192E-1	8.331E-1	25.85	29.01	22.90	22.96	21.69	25.90	14.41	26.94	6.07	19.03	15.78	17.24	40.70	14.43	0.14	21.57	26.48	34.02	22.31	25.93	5.243E-5	2.667E-5	1.061E-4	1.044E-4	1.452E-4	5.178E-5	1.473E-3	4.118E-5	9.752E-2	3.079E-4	8.873E-4	5.41E-4	3.607E-6	1.458E-3	10E-1	1.501E-4	4.551E-5	1.042E-5	1.233E-4	5.147E-5	KR	ACDEFGHILMNPQSTVWY
5	12.92	2.672E-3	'TTTTTTTTTNSN	0.0	N	9.99E-1	10E-1	9.999E-1	6.877E-1	5E-1	10E-1	17.18	4.88	24.32	28.99	12.32	9.13	12.69	13.14	26.40	15.82	16.45	2.00	33.07	10.41	27.53	2.52	0.84	8.53	15.33	9.30	5.518E-4	2.105E-1	7.478E-5	2.677E-5	3.446E-3	1.609E-2	2.938E-3	2.438E-3	4.631E-5	8.773E-4	7.03E-4	8.952E-1	1.233E-5	8.339E-3	3.626E-5	7.682E-1	9.984E-1	2.227E-2	1.045E-3	1.471E-2	CGNSTVY	ADEFHIKLMPQRW
6	11.47	5.031E-3	'RRHHH-RHR---	0.473105378924215	N	3.057E-1	9.785E-1	6.921E-1	8.888E-1	9.994E-1	9.999E-1	15.74	16.60	12.16	9.21	15.29	12.00	2.73	19.98	4.65	16.59	11.16	9.18	20.33	5.16	2.46	8.84	10.97	21.74	9.77	11.78	9.028E-4	6.695E-4	3.691E-3	1.544E-2	1.058E-3	3.966E-3	7.079E-1	2.329E-4	2.441E-1	6.717E-4	5.816E-3	1.569E-2	2.107E-4	1.754E-1	7.847E-1	1.882E-2	6.351E-3	1.432E-4	1.147E-2	4.381E-3	EHKNQRSW	ACDFGILMPTVY
"""


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
            'monospace', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL);
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
            'monospace', cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL);
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

def get_response(fs):
    """
    @param fs: a FieldStorage object containing the cgi arguments
    @return: a (response_headers, response_text) pair
    """
    # start writing the response type
    response_headers = []
    # read the headers
    headers = Util.get_stripped_lines(StringIO(fs.headers))
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
