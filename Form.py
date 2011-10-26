"""
Mediate between the interface and the scripts.

Note that at least three interfaces are involved.
The first is a very simple web-based interface which was originally
based on a simple cgi subset of wsgi and has been replaced
by a very simple cherrypy web-based interface.
The second interface is through the galaxy project.
The third interface is through the mobyle project.
"""

import cgi
from StringIO import StringIO
import itertools
import string

import numpy as np

import MatrixUtil
import iterutils
import Util

try:
    from lxml import etree
except ImportError as e:
    pass

g_safe_letters = string.letters + string.digits + '_.-'

g_imageformat_to_contenttype = {
        'svg' : 'image/svg+xml',
        'png' : 'image/png',
        'pdf' : 'application/pdf',
        'ps' : 'application/postscript'}

g_imageformat_to_ext = {
        'svg' : 'svg',
        'png' : 'png',
        'pdf' : 'pdf',
        'ps' : 'ps'}

g_imageformat_to_r_function = {
        'svg' : 'svg',
        'png' : 'png',
        'pdf' : 'pdf',
        'ps' : 'postscript'}

g_tikzformat_to_contenttype = {
        'tikz' : 'text/plain',
        'tex' : 'text/plain',
        'pdf' : 'application/pdf',
        'png' : 'image/png'}

g_latexformat_to_contenttype = {
        'tex' : 'text/plain',
        'pdf' : 'application/pdf',
        'png' : 'image/png'}

g_tikzformat_to_ext = {
        'tikz' : 'tikz',
        'tex' : 'tex',
        'pdf' : 'pdf',
        'png' : 'png'}

g_latexformat_to_ext = {
        'tex' : 'tex',
        'pdf' : 'pdf',
        'png' : 'png'}

class HelpItem(object):
    def __init__(self, command, description):
        self.command = command
        self.description = description
    def is_isolated(self):
        return False
    def _get_partial_line(self, depth):
        return (' ' * 2*depth) + self.command
    def _get_partial_lines(self, depth):
        return [self._get_partial_line(depth)]
    def get_partial_lines(self):
        return self._get_partial_lines(0)
    def get_full_lines(self, depth, max_partial_length):
        elements = [
                self._get_partial_line(depth).ljust(max_partial_length),
                '  : ',
                self.description]
        return [''.join(elements)]

class HelpGroup(object):
    def __init__(self, description):
        self.isolated = True
        self.description = description
        self.subitems = []
    def is_isolated(self):
        return True
    def _get_partial_lines(self, depth):
        lines = []
        for item in self.subitems:
            lines.extend(item._get_partial_lines(depth+1))
        return lines
    def get_partial_lines(self):
        return self._get_partial_lines(0)
    def get_full_lines(self, depth, max_partial_length):
        lines = [(' ' * 2*depth) + self.description + ':']
        for item in self.subitems:
            lines.extend(item.get_full_lines(depth+1, max_partial_length))
        return lines

def any_isolated(seq):
    return any(x.is_isolated() for x in seq)

def get_help_string(form_objects):
    """
    @param form_objects: a list of form objects
    """
    # Get the help objects
    help_objects = [
            x.get_help_object() for x in form_objects if not x.web_only()]
    # Get the length of the longest partial line.
    pline_lists = [x.get_partial_lines() for x in help_objects]
    plines = list(itertools.chain.from_iterable(pline_lists))
    max_len = Util.max_length(plines)
    # True if a blank line should be added after the corresponding help object.
    bsep = [any_isolated(p) for p in iterutils.pairwise(help_objects)] + [0]
    # Add these lines after the corresponding help object.
    lsep = [[''] if x else [] for x in bsep]
    # Build the list of output lines.
    lines = []
    for obj, sep in zip(help_objects, lsep):
        lines.extend(obj.get_full_lines(0, max_len))
        lines.extend(sep)
    return '\n'.join(lines)

def get_html_string(form_objects):
    """
    @param form_objects: a list of form objects
    """
    lines = []
    for form_object in form_objects:
        lines.extend(form_object.get_html_lines())
        if form_object is not form_objects[-1]:
            lines.append('<br/><br/>')
    return '\n'.join(lines)


class FormError(Exception):
    pass


def _get_checkbox_line(esc_label, checked):
    """
    @param esc_label: escaped label
    @param checked: boolean to specify whether item is checked or not
    """
    lines = (
            'input type="checkbox" name="%s"' % esc_label,
            'id="%s"' % esc_label,
            'value="%s"' % esc_label)
    base = ' '.join(lines)
    if checked:
        return '<%s checked="yes"/>' % base
    else:
        return '<%s/>' % base

def _get_radio_line(esc_group_label, esc_label, checked):
    """
    @param esc_group_label: escaped label for the whole group
    @param esc_label: escaped label for the particular button
    @param checked: boolean to specify whether item is checked or not
    """
    lines = (
            'input type="radio" name="%s"' % esc_group_label,
            'id="%s"' % esc_label,
            'value="%s"' % esc_label)
    base = ' '.join(lines)
    if checked:
        return '<%s checked="yes"/>' % base
    else:
        return '<%s/>' % base

def _get_colon_label_line(esc_label, esc_description):
    """
    @param esc_label: escaped label
    @param esc_description: escaped description
    @return: labeled description
    """
    return '<label for="%s">%s:</label>' % (esc_label, esc_description)

def _get_label_line(esc_label, esc_description):
    """
    @param esc_label: escaped label
    @param esc_description: escaped description
    @return: labeled description
    """
    return '<label for="%s">%s</label>' % (esc_label, esc_description)

def _get_textbox_line(esc_label, esc_default_line, width):
    """
    @param esc_label: escaped label
    @param esc_default_line: escaped default line
    @param width: width of the textbox
    """
    lines = (
            '<input type="text" name="%s"' % esc_label,
            'id="%s"' % esc_label,
            'value="%s"' % esc_default_line,
            'size="%d"/>' % width)
    return ' '.join(lines)

def _get_textarea_header(esc_label, nrows):
    """
    @param esc_label: escaped label
    @param nrows: the number of rows in the text area
    @return: the single line header for the textarea
    """
    lines = (
            '<textarea name="%s"' % esc_label,
            'id="%s"' % esc_label,
            'rows="%d" cols="70" wrap="off">' % nrows)
    return ' '.join(lines)

def _set_unique(d_out, label, value):
    if label in d_out:
        raise FormError(label + ' is duplicated')
    d_out[label] = value

class RadioGroup:
    """
    In html, radio buttons must be grouped.
    """

    def __init__(self, label, description, radio_items):
        """
        @param label: something like a variable name
        @param description: single line description of the group
        @param radio_items: a sequence of RadioItem objects
        """
        # initialize the member variables
        self.label = label
        self.description = description
        self.radio_items = radio_items
        # assert that exactly one radio_item is checked
        self._get_checked_list()

    def web_only(self):
        return False

    def _get_checked_list(self):
        checked_list = [item for item in self.radio_items if item.default]
        if len(checked_list) != 1:
            msg = 'exactly one radio button should be checked by default'
            raise FormError(msg)
        return checked_list
    
    def _get_default_string(self):
        return self._get_checked_list()[0].label

    def get_galaxy_cmd(self):
        return '--%s=$%s' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        param = etree.SubElement(parent, 'param',
                type='select', display='radio', multiple='false',
                name=self.label, label=self.description)
        for item in self.radio_items:
            if item.default:
                etree.SubElement(param, 'option',
                        value=item.label,
                        selected='true').text = item.description
            else:
                etree.SubElement(param, 'option',
                        value=item.label).text = item.description

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        mobyle_class = 'Choice'
        meta_code = '" --" + str(value)'
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        vdef = etree.SubElement(parameter, 'vdef')
        etree.SubElement(vdef, 'value').text = self._get_default_string()
        vlist = etree.SubElement(parameter, 'vlist')
        for item in self.radio_items:
            velem = etree.SubElement(vlist, 'velem')
            etree.SubElement(velem, 'value').text = item.label
            etree.SubElement(velem, 'label').text = item.description
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        default_string = self._get_default_string()
        description = '%s (default: %s)' % (self.description, default_string)
        obj = HelpGroup(description)
        obj.subitems = [x.get_help_object() for x in self.radio_items]
        return obj

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        escaped_description = cgi.escape(self.description)
        # add the line that describes the group
        description_line = escaped_description + ':'
        lines.append(description_line)
        lines.append('<br/>')
        # add the lines that describe the entries
        for radio_item in self.radio_items:
            lines.extend(radio_item.get_html_lines(self.label))
            if radio_item is not self.radio_items[-1]:
                lines.append('<br/>')
        # return the list of lines
        return lines

    def process_cmdline_dict(self, d_in, d_out):
        selected_item = None
        for item in self.radio_items:
            if item.label in d_in:
                if d_in[item.label] != True:
                    msg_a = 'to select the %s option ' % item.label
                    msg_b = 'use --%s' % item.label
                    raise FormError(msg_a + msg_b)
                if selected_item:
                    msg_a = 'multiple radio button selections: '
                    msg_b = '%s and %s' % (selected_item.label, item.label)
                    raise FormError(msg_a + msg_b)
                selected_item = item
        if selected_item is None:
            for item in self.radio_items:
                if item.default:
                    if selected_item:
                        msg = 'multiple radio button default selections'
                        raise FormError(msg)
                    selected_item = item
        if selected_item is None:
            msg = 'no default radio button selection'
            raise FormError(msg)
        for item in self.radio_items:
            _set_unique(d_out, item.label, item is selected_item)
        _set_unique(d_out, self.label, selected_item.label)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            lines = (
                    'no radio button option',
                    'was selected for the field "%s"' % self.label)
            raise FormError(' '.join(lines))
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        # Assert that the selected value
        # is actually one of the radio button options.
        if value not in set(item.label for item in self.radio_items):
            msg = 'an invalid radio button option was selected: %s' % value
            raise FormError(msg)
        # For the group,
        # set the value for the attribute in the fieldstorage object.
        setattr(fs, self.label, value)
        # For each item,
        # set the value for the attribute in the fieldstorage object.
        for radio_item in self.radio_items:
            # verify that the attribute is not already taken
            if hasattr(fs, radio_item.label):
                msg = 'the object already has the attribute "%s"' % self.label
                raise FormError(msg)
            is_checked = (radio_item.label == value)
            setattr(fs, radio_item.label, is_checked)

class CheckGroup:
    """
    In html, checkboxes are not grouped together.
    This extra level of abstraction is to allow
    sets of checkboxes to be visually separated from each other
    and from the rest of the page.
    """

    def __init__(self, label, description, check_items):
        """
        @param label: something like a variable name
        @param description: single line description of the group
        @param check_items: a sequence of CheckItem objects
        """
        self.label = label
        self.description = description
        self.check_items = check_items

    def web_only(self):
        return False

    def get_galaxy_cmd(self):
        return '--%s="$%s"' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        param = etree.SubElement(parent, 'param',
                type='select', display='checkboxes', multiple='true',
                name=self.label, label=self.description)
        for item in self.check_items:
            if item.default:
                etree.SubElement(param, 'option',
                        value=item.label,
                        selected='true').text = item.description
            else:
                etree.SubElement(param, 'option',
                        value=item.label).text = item.description

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        paragraph = etree.SubElement(parent, 'paragraph')
        etree.SubElement(paragraph, 'name').text = self.label
        etree.SubElement(paragraph, 'prompt').text = self.description
        parameters = etree.SubElement(paragraph, 'parameters')
        for item in self.check_items:
            item.add_mob_xml(parameters, next_argpos)
            next_argpos += 1
        return len(self.check_items)

    def get_help_object(self):
        obj = HelpGroup(self.description)
        obj.subitems = [x.get_help_object() for x in self.check_items]
        return obj

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        escaped_description = cgi.escape(self.description)
        # add the line that describes the group
        description_line = escaped_description + ':'
        lines.append(description_line)
        lines.append('<br/>')
        # add the lines that describe the entries
        for check_item in self.check_items:
            lines.extend(check_item.get_html_lines())
            if check_item is not self.check_items[-1]:
                lines.append('<br/>')
        # return the list of lines
        return lines

    def process_cmdline_dict(self, d_in, d_out):
        for item in self.check_items:
            item.process_cmdline_dict(d_in, d_out)
        checked = set(x.label for x in self.check_items if x.label in d_out)
        # look for items set using the checkgroup itself under galaxy
        if self.label in d_in:
            observed = set(d_in[self.label].split(','))
            allowed = set(x.label for x in self.check_items)
            bad = observed - allowed
            if bad:
                raise FormError('these choices are not allowed: %s' % bad)
            for item in self.check_items:
                d_out[item.label] = True
            checked |= observed
        _set_unique(d_out, self.label, checked)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # Decorate the FieldStorage object
        # with the boolean values of each item separately.
        for check_item in self.check_items:
            check_item.process_fieldstorage(fs)
        # get the set of checked item labels
        checked_set = set()
        for check_item in self.check_items:
            values = fs.getlist(check_item.label)
            if values:
                checked_set.add(check_item.label)
        # decorate the FieldStorage object with the set of checked selections
        setattr(fs, self.label, checked_set)


class RadioItem:

    def __init__(self, label, description, default=False):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default: True if the group item is selected
        """
        self.label = label
        self.description = description
        self.default = default

    def web_only(self):
        return False

    def get_help_object(self):
        return HelpItem(
                '--%s' % self.label,
                self.description)

    def get_html_lines(self, group_label):
        """
        @param group_label: the label of the radio button group
        @return: the list of lines of html text
        """
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_group_label = cgi.escape(group_label)
        lines = [
                _get_radio_line(esc_group_label, esc_label, self.default),
                _get_label_line(esc_label, esc_description)]
        return lines


class CheckItem:

    def __init__(self, label, description, default=False):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default: True if the group item is selected
        """
        self.label = label
        self.description = description
        self.default = default

    def web_only(self):
        return False

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        vdef_text = '1' if self.default else '0'
        mobyle_class = 'Boolean'
        meta_code = '" --%s=" + ("Y" if value else "N")' % self.label
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        vdef = etree.SubElement(parameter, 'vdef')
        etree.SubElement(vdef, 'value').text = vdef_text
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        s = 'Y' if self.default else 'N'
        return HelpItem(
                '--%s=%s' % (self.label, s),
                self.description)

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        escaped_label = cgi.escape(self.label)
        escaped_description = cgi.escape(self.description)
        lines = [
                _get_checkbox_line(escaped_label, self.default),
                _get_label_line(escaped_label, escaped_description)]
        return lines

    def process_cmdline_dict(self, d_in, d_out):
        value_in = d_in.get(self.label, None)
        if value_in is None:
            value_out = self.default
        elif value_in in (True, 'Y'):
            value_out = True
        elif value_in == 'N':
            value_out = False
        else:
            raise FormError('bad checkbox value: ' + value_in)
        _set_unique(d_out, self.label, value_out)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        # set the value for the attribute in the fieldstorage object
        is_checked = True if values else False
        setattr(fs, self.label, is_checked)


class SingleLine:
    """
    This represents a single line text box.
    """

    def __init__(self, label, description, default_line):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_line: the default line of text
        """
        self.label = label
        self.description = description
        self.default_line = default_line

    def web_only(self):
        return False

    def get_galaxy_cmd(self):
        return '--%s="$%s"' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        etree.SubElement(parent, 'param',
                type='text', value=self.default_line,
                name=self.label, label=self.description)

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        mobyle_class = 'String'
        vdef_text = cgi.escape(self.default_line)
        meta_code = '" --%s=" + str(value)' % self.label
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        vdef = etree.SubElement(parameter, 'vdef')
        etree.SubElement(vdef, 'value').text = vdef_text
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        if set(self.default_line) - set(g_safe_letters):
            s = '"%s"' % self.default_line
        else:
            s = self.default_line
        return HelpItem(
                '--%s=%s' % (self.label, s),
                self.description)

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        # calculate a multiple of ten that will hold the string
        width = ((len(self.default_line) / 10) + 1) * 10
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_line = cgi.escape(self.default_line)
        lines = [
                # add the label line followed by a line break
                _get_colon_label_line(esc_label, esc_description),
                '<br/>',
                # add the textbox line
                _get_textbox_line(esc_label, esc_default_line, width)]
        # return the list of lines
        return lines

    def process_cmdline_dict(self, d_in, d_out):
        value = d_in.get(self.label, self.default_line)
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Float:
    """
    An floating point number is requested using a single line in the form.
    """

    def __init__(self, label, description, default_float,
            low_exclusive=None, low_inclusive=None,
            high_exclusive=None, high_inclusive=None):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_float: the default floating point number
        """
        # assert that at most one low bound and at most one high bound is set
        if (low_exclusive is not None) and (low_inclusive is not None):
            raise ValueError('both low_exclusive and low_inclusive were set')
        if (high_exclusive is not None) and (high_inclusive is not None):
            raise ValueError('both high_exclusive and high_inclusive were set')
        # set the member variables
        self.label = label
        self.description = description
        self.default_float = default_float
        self.low_exclusive = low_exclusive
        self.low_inclusive = low_inclusive
        self.high_exclusive = high_exclusive
        self.high_inclusive = high_inclusive

    def web_only(self):
        return False

    def get_galaxy_cmd(self):
        return '--%s=$%s' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        etree.SubElement(parent, 'param',
                type='float', value=str(self.default_float),
                name=self.label, label=self.description)

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        mobyle_class = 'Float'
        vdef_text = str(self.default_float)
        meta_code = '" --%s=" + str(value)' % self.label
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        vdef = etree.SubElement(parameter, 'vdef')
        etree.SubElement(vdef, 'value').text = vdef_text
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        return HelpItem(
                '--' + self.label + '=' + str(self.default_float),
                self.description)

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        # get the default line of text
        default_line = str(self.default_float)
        # calculate a multiple of ten that will hold the string
        width = ((len(default_line) / 10) + 1) * 10
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_line = cgi.escape(default_line)
        lines = [
                # add the label line followed by a line break
                _get_colon_label_line(esc_label, esc_description),
                '<br/>',
                # add the textbox line
                _get_textbox_line(esc_label, esc_default_line, width)]
        # return the list of lines
        return lines

    def validate(self, value):
        identifier = 'the floating point number in the field "%s"' % self.label
        if self.low_exclusive is not None:
            if value <= self.low_exclusive:
                lines = (
                        '%s must be' % identifier,
                        'greater than %f' % self.low_exclusive)
                raise FormError(' '.join(lines))
        if self.low_inclusive is not None:
            if value < self.low_inclusive:
                lines = (
                        '%s must be' % identifier,
                        'greater than or equal to %f' % self.low_inclusive)
                raise FormError(' '.join(lines))
        if self.high_exclusive is not None:
            if value >= self.high_exclusive:
                lines = (
                        '%s must be' % identifier,
                        'less than %f' % self.high_exclusive)
                raise FormError(' '.join(lines))
        if self.high_inclusive is not None:
            if value > self.high_inclusive:
                lines = (
                        '%s must be' % identifier,
                        'less than or equal to %f' % self.high_inclusive)
                raise FormError(' '.join(lines))

    def process_cmdline_dict(self, d_in, d_out):
        value_s = d_in.get(self.label, None)
        if value_s is None:
            value = float(self.default_float)
        else:
            try:
                value = float(value_s)
            except ValueError as e:
                raise FormError(value_s + ' is not a floating point value')
        self.validate(value)
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            lines = (
                    'no floating point number',
                    'was given for the field "%s"' % self.label)
            raise FormError(' '.join(lines))
        elif len(values) == 1:
            value_string = values[0]
            try:
                value = float(value_string)
            except ValueError, e:
                lines = (
                        '%s could not be interpreted' % value_string,
                        'as a floating point number')
                raise FormError(' '.join(lines))
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        self.validate(value)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Integer:
    """
    An integer is requested using a single line in the form.
    """

    def __init__(self, label, description, default_integer,
            low=None, high=None):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_integer: the default integer
        @param low: the integer must be greater than or equal to this value
        @param high: the integer must be less than or equal to this value
        """
        self.label = label
        self.description = description
        self.default_integer = default_integer
        self.low = low
        self.high = high

    def web_only(self):
        return False

    def get_galaxy_cmd(self):
        return '--%s=$%s' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        etree.SubElement(parent, 'param',
                type='integer', value=str(self.default_integer),
                name=self.label, label=self.description)

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        mobyle_class = 'Integer'
        vdef_text = str(self.default_integer)
        meta_code = '" --%s=" + str(value)' % self.label
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        vdef = etree.SubElement(parameter, 'vdef')
        etree.SubElement(vdef, 'value').text = vdef_text
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        return HelpItem(
                '--' + self.label + '=' + str(self.default_integer),
                self.description)

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        # get the default line of text
        default_line = str(self.default_integer)
        # calculate a multiple of ten that will hold the string
        width = ((len(default_line) / 10) + 1) * 10
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_line = cgi.escape(default_line)
        lines = [
                # add the label line followed by a line break
                _get_colon_label_line(esc_label, esc_description),
                '<br/>',
                # add the textbox line
                _get_textbox_line(esc_label, esc_default_line, width)]
        # return the list of lines
        return lines

    def validate(self, value):
        if self.low is not None:
            if value < self.low:
                msg_a = 'the integer in the field "%s" ' % self.label
                msg_b = 'must be at least %d' % self.low
                raise FormError(msg_a + msg_b)
        if self.high is not None:
            if value > self.high:
                msg_a = 'the integer in the field "%s" ' % self.label
                msg_b = 'must be at most %d' % self.high
                raise FormError(msg_a + msg_b)

    def process_cmdline_dict(self, d_in, d_out):
        value_s = d_in.get(self.label, None)
        if value_s is None:
            value = int(self.default_integer)
        else:
            try:
                value = int(value_s)
            except ValueError as e:
                raise FormError(value_s + ' is not an integer')
        self.validate(value)
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            msg = 'no integer was given for the field "%s"' % self.label
            raise FormError(msg)
        elif len(values) == 1:
            value_string = values[0]
            try:
                value = int(value_string)
            except ValueError, e:
                raise FormError('%s is not an integer' % value_string)
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        self.validate(value)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Matrix:
    """
    This represents a numpy matrix entered by the user.
    The matrix should be entered in a multi-line text box.
    Each non-empty line should be a row of whitespace separated numbers.
    """

    def __init__(self, label, description, default_matrix,
            matrix_assertion=None):
        """
        The matrix_assertion may be a function or None.
        If it is a function then it may raise MatrixUtil.MatrixError.
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_matrix: the default numpy array
        @param matrix_assertion: None or a function of a matrix
        """
        self.label = label
        self.description = description
        self.default_matrix = default_matrix
        self.matrix_assertion = matrix_assertion

    def web_only(self):
        return False

    def get_galaxy_cmd(self):
        return '--%s=$%s' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        etree.SubElement(parent, 'param',
                type='data', format='txt',
                name=self.label, label=self.description)

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        mobyle_class = 'Text'
        matrix_string = MatrixUtil.m_to_string(self.default_matrix)
        example_text = cgi.escape(matrix_string)
        meta_code = '" --%s=" + str(value)' % self.label
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        example = etree.SubElement(parameter, 'example').text = example_text
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        return HelpItem(
                '--' + self.label + '=FILENAME',
                self.description)

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        # get the number of rows to use for the textarea
        sio = StringIO(MatrixUtil.m_to_string(self.default_matrix))
        nrows = len(list(sio.readlines())) + 1
        nrows = min(nrows, 12)
        # get the matrix as an unescaped string
        default_string = MatrixUtil.m_to_string(self.default_matrix)
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_string = cgi.escape(default_string)
        lines = [
                # add the label line followed by a line break
                _get_colon_label_line(esc_label, esc_description),
                '<br/>',
                # add the textarea
                _get_textarea_header(esc_label, nrows),
                esc_default_string,
                '</textarea>']
        # return the list of lines
        return lines

    def process_cmdline_dict(self, d_in, args_out):
        filename = d_in.get(self.label, None)
        if filename is None:
            value = self.default_matrix
        else:
            with open(filename) as fin:
                try:
                    value = np.array(MatrixUtil.read_matrix(fin))
                    if self.matrix_assertion:
                        self.matrix_assertion(value)
                except MatrixUtil.MatrixError, e:
                    raise FormError(e)
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            msg = 'the value for the field "%s" is empty' % self.label
            raise FormError(msg)
        elif len(values) == 1:
            try:
                value = np.array(MatrixUtil.read_matrix(StringIO(values[0])))
                if self.matrix_assertion:
                    self.matrix_assertion(value)
            except MatrixUtil.MatrixError, e:
                raise FormError(e)
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class MultiLine:
    """
    This represents a multi-line text box.
    The text wrapping is always off,
    and the number of columns is always seventy.
    """

    def __init__(self, label, description, default_string):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_string: the default multi-line string of text
        """
        self.label = label
        self.description = description
        self.default_string = default_string

    def get_default_string(self):
        return self.default_string

    def web_only(self):
        return False

    def get_galaxy_cmd(self):
        return '--%s=$%s' % (self.label, self.label)

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        etree.SubElement(parent, 'param',
                type='data', format='txt',
                name=self.label, label=self.description)

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        mobyle_class = 'Text'
        example_text = cgi.escape(self.get_default_string())
        meta_code = '" --%s=" + str(value)' % self.label
        parameter = etree.SubElement(
                parent, 'parameter', ismandatory='1', issimple='1')
        etree.SubElement(parameter, 'name').text = self.label
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = self.description
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        etree.SubElement(datatype, 'class').text = mobyle_class
        fmt = etree.SubElement(parameter, 'format')
        code = etree.SubElement(fmt, 'code', proglang='python')
        code.text = meta_code
        example = etree.SubElement(parameter, 'example').text = example_text
        etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
        return 1

    def get_help_object(self):
        return HelpItem(
                '--' + self.label + '=FILENAME',
                self.description)

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        # get the number of rows to use for the textarea
        sio = StringIO(self.get_default_string())
        nrows = len(list(sio.readlines())) + 1
        nrows = min(nrows, 12)
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_string = cgi.escape(self.get_default_string())
        lines = [
                # add the label line followed by a line break
                _get_colon_label_line(esc_label, esc_description),
                '<br/>',
                # add the textarea
                _get_textarea_header(esc_label, nrows),
                esc_default_string,
                '</textarea>']
        # return the list of lines
        return lines

    def process_cmdline_dict(self, d_in, d_out):
        filename = d_in.get(self.label, None)
        if filename is None:
            value = self.get_default_string()
        else:
            with open(filename) as fin:
                value = fin.read()
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Sequence(MultiLine):
    """
    This represents a tuple of strings.
    It is not supposed to be a DNA sequence or something.
    Whitespace at either end of a string is ignored.
    Strings consisting entirely of whitespace are ignored.
    """

    def __init__(self, label, description, default_list):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_list: the default list of strings
        """
        self.label = label
        self.description = description
        self.default_list = tuple(self._gen_reduced(default_list))

    def _gen_reduced(self, arr):
        for s in arr:
            rs = s.strip()
            if rs:
                yield rs

    def get_default_string(self):
        return '\n'.join(self.default_list)

    def process_cmdline_dict(self, d_in, d_out):
        filename = d_in.get(self.label, None)
        if filename is None:
            value = self.default_list
        else:
            with open(filename) as fin:
                value = tuple(self._gen_reduced(fin.readlines()))
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            msg = 'the object already has the attribute "%s"' % self.label
            raise FormError(msg)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            msg = 'the value for the field "%s" is ambiguous' % self.label
            raise FormError(msg)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, tuple(self._gen_reduced(value.splitlines())))


class ContentDisposition(RadioGroup):
    def __init__(self):
        """
        The group label is contentdisposition.
        The description and radio items are hard coded.
        """
        label = 'contentdisposition'
        description = 'delivery options'
        radio_items = [
                RadioItem('inline', 'view', True),
                RadioItem('attachment', 'download')]
        RadioGroup.__init__(self, label, description, radio_items)

    def web_only(self):
        return True


class ImageFormat(RadioGroup):
    def __init__(self):
        """
        The group label is imageformat.
        The description and radio items are hard coded.
        """
        label = 'imageformat'
        description = 'output image format'
        radio_items = [
            RadioItem('png', 'png', True),
            RadioItem('pdf', 'pdf'),
            RadioItem('ps', 'postscript'),
            RadioItem('svg', 'svg')]
        RadioGroup.__init__(self, label, description, radio_items)


class TikzFormat(RadioGroup):
    def __init__(self):
        """
        The group label is tikzformat.
        The description and radio items are hard coded.
        """
        label = 'tikzformat'
        description = 'output format'
        radio_items = [
            RadioItem('tikz', 'TikZ code', True),
            RadioItem('tex', 'full LaTeX code'),
            RadioItem('pdf', 'pdf'),
            RadioItem('png', 'png')]
        RadioGroup.__init__(self, label, description, radio_items)


class LatexFormat(RadioGroup):
    def __init__(self):
        """
        The group label is latexformat.
        The description and radio items are hard coded.
        """
        label = 'latexformat'
        description = 'output format'
        radio_items = [
            RadioItem('tex', 'LaTeX code', True),
            RadioItem('pdf', 'pdf'),
            RadioItem('png', 'png')]
        RadioGroup.__init__(self, label, description, radio_items)

