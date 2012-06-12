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

import smallutil
import matrixio
from matrixio import MatrixIOError

try:
    from lxml import etree
except ImportError as e:
    pass

# reduced from 70
g_default_textarea_ncols = 50

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

class Preset(object):
    def __init__(self, description, d):
        self.description = description
        self.d = d

def get_default_preset(form_objects):
    preset_pairs = []
    for obj in form_objects:
        obj_preset_pairs = None
        try:
            obj_preset_pairs = obj.get_preset_pairs()
        except AttributeError as e:
            pass
        if obj_preset_pairs is None:
            obj_preset_pairs = [obj.get_preset_pair()]
        preset_pairs.extend(obj_preset_pairs)
    return Preset('default', dict(preset_pairs))

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
    help_objects = []
    for obj in form_objects:
        obj_help_objects = None
        try:
            obj_help_objects = obj.get_help_objects()
        except AttributeError as e:
            pass
        if obj_help_objects is None:
            obj_help_objects = [obj.get_help_object()]
        help_objects.extend(obj_help_objects)
    # Get the length of the longest partial line.
    pline_lists = [x.get_partial_lines() for x in help_objects]
    plines = list(itertools.chain.from_iterable(pline_lists))
    max_len = max(len(x) for x in plines)
    # True if a blank line should be added after the corresponding help object.
    bsep = [any_isolated(p) for p in smallutil.pairwise(help_objects)] + [0]
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


#TODO include group label
def _get_checkbox_line(esc_label, checked):
    """
    @param esc_label: escaped label
    @param checked: boolean to specify whether item is checked or not
    """
    lines = (
            'input type="checkbox"',
            'name="%s"' % esc_label,
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
            'rows="%d"' % nrows,
            'cols="%d"' % g_default_textarea_ncols,
            'wrap="off">')
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

    def get_preset_pair(self):
        return self.label, self._get_checked_list()[0].label

    def _get_checked_list(self):
        checked_list = [item for item in self.radio_items if item.default]
        if len(checked_list) != 1:
            raise FormError(
                    'exactly one radio button should be checked by default')
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
                    raise FormError(
                            'to select the %s option '
                            'use --%s' % (item.label, item.label))
                if selected_item:
                    raise FormError(
                            'multiple radio button selections: '
                            '%s and %s' % (selected_item.label, item.label))
                selected_item = item
        if selected_item is None:
            for item in self.radio_items:
                if item.default:
                    if selected_item:
                        raise FormError(
                                'multiple radio button default selections')
                    selected_item = item
        if selected_item is None:
            raise FormError(
                    'no default radio button selection')
        for item in self.radio_items:
            _set_unique(d_out, item.label, item is selected_item)
        _set_unique(d_out, self.label, selected_item.label)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError(
                    'no radio button option '
                    'was selected for the field "%s"' % self.label)
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
        # Assert that the selected value
        # is actually one of the radio button options.
        if value not in set(item.label for item in self.radio_items):
            raise FormError(
                    'an invalid radio button option was selected: %s' % value)
        # For the group,
        # set the value for the attribute in the fieldstorage object.
        setattr(fs, self.label, value)
        # For each item,
        # set the value for the attribute in the fieldstorage object.
        for radio_item in self.radio_items:
            # verify that the attribute is not already taken
            if hasattr(fs, radio_item.label):
                raise FormError(
                        'the object already has '
                        'the attribute "%s"' % self.label)
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

    def get_preset_pair(self):
        values = tuple(item.label for item in self.check_items if item.default)
        return self.label, values

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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
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

    def get_preset_pair(self):
        return self.label, self.default_line

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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
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

    def get_preset_pair(self):
        return self.label, self.default_float

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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError(
                    'no floating point number '
                    'was given for the field "%s"' % self.label)
        elif len(values) == 1:
            value_string = values[0]
            try:
                value = float(value_string)
            except ValueError as e:
                raise FormError(
                        '%s could not be interpreted '
                        'as a floating point number' % value_string)
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
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

    def get_preset_pair(self):
        return self.label, self.default_integer

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
                raise FormError(
                        'the integer in the field "%s" '
                        'must be at least %d' % (self.label, self.low))
        if self.high is not None:
            if value > self.high:
                raise FormError(
                        'the integer in the field "%s" '
                        'must be at most %d' % (self.label, self.high))

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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError(
                    'no integer was given for the field "%s"' % self.label)
        elif len(values) == 1:
            value_string = values[0]
            try:
                value = int(value_string)
            except ValueError as e:
                raise FormError('%s is not an integer' % value_string)
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
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
        If it is a function then it may raise matrixio.MatrixIOError.
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

    def get_preset_pair(self):
        return self.label, self.default_matrix

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
        matrix_string = matrixio.m_to_string(self.default_matrix)
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
        sio = StringIO(matrixio.m_to_string(self.default_matrix))
        nrows = len(list(sio.readlines())) + 1
        nrows = min(nrows, 12)
        # get the matrix as an unescaped string
        default_string = matrixio.m_to_string(self.default_matrix)
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

    def process_cmdline_dict(self, d_in, d_out):
        #FIXME this had a typo and has apparently never been used successfully
        filename = d_in.get(self.label, None)
        if filename is None:
            value = self.default_matrix
        else:
            with open(filename) as fin:
                try:
                    value = np.array(matrixio.read_matrix(fin))
                    if self.matrix_assertion:
                        self.matrix_assertion(value)
                except MatrixIOError as e:
                    raise FormError(e)
        _set_unique(d_out, self.label, value)

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError(
                    'the value for the field "%s" is empty' % self.label)
        elif len(values) == 1:
            try:
                value = np.array(matrixio.read_matrix(StringIO(values[0])))
                if self.matrix_assertion:
                    self.matrix_assertion(value)
            except MatrixIOError as e:
                raise FormError(e)
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class MultiLine:
    """
    This represents a multi-line text box.
    The text wrapping is always off, and the number of columns is fixed.
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

    def get_preset_pair(self):
        return self.label, self.get_default_string()

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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)

class _Interval:

    def web_only(self):
        return False


class IntegerInterval(_Interval):
    """
    This represents an ordered pair of integers.
    It is experimental.
    """

    def __init__(self,
            first_label, second_label,
            description,
            first_default, second_default,
            low=None, high=None, low_width=None, high_width=None):
        """
        @param first_label: something like a variable name
        @param second_label: something like a variable name
        @param description: a single line description of the item
        @param first_default: the default low integer
        @param second_default: the default high integer
        @param low: inclusive lower bound on the first and second values
        @param high: inclusive upper bound on the first and second values
        @param low_width: inclusive lower bound on width
        @param high_width: inclusive upper bound on width
        """
        self.first_label = first_label
        self.second_label = second_label
        self.description = description
        self.first_default = first_default
        self.second_default = second_default
        self.low = low
        self.high = high
        self.low_width = low_width
        self.high_width = high_width

    def _get_temp_items(self):
        """
        Make a couple of temporary Integer items for convenience.
        """
        return (
                Integer(self.first_label, self.description + ' (low)',
                    self.first_default, low=self.low, high=self.high),
                Integer(self.second_label, self.description + ' (high)',
                    self.second_default, low=self.low, high=self.high),
                )

    def get_preset_pairs(self):
        return tuple(item.get_preset_pair() for item in self._get_temp_items())

    def get_galaxy_cmd(self):
        a, b = self._get_temp_items()
        return '%s %s' % (a.get_galaxy_cmd(), b.get_galaxy_cmd())

    def add_galaxy_xml(self, parent):
        """
        Add a galaxy parameter to the xml tree.
        @param parent: parent etree element
        """
        for item in self._get_temp_items():
            item.add_galaxy_xml(parent)

    def add_mob_xml(self, parent, next_argpos):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param next_argpos: a 1-based integer for cmdline arg ordering
        @return: the number of args added on the command line
        """
        items = self._get_temp_items()
        argpos = next_argpos
        for i, item in items:
            argpos = item.add_mob_xml(parent, argpos)
        return len(items)

    def get_help_objects(self):
        return tuple(item.get_help_object() for item in self._get_temp_items())

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        # get the default line of text
        default_line_a = str(self.first_default)
        default_line_b = str(self.second_default)
        # calculate a multiple of ten that will hold the string
        width = max(
                ((len(default_line_a) / 10) + 1) * 10,
                ((len(default_line_b) / 10) + 1) * 10,
                )
        # get escaped values
        esc_label_a = cgi.escape(self.first_label)
        esc_label_b = cgi.escape(self.second_label)
        esc_description = cgi.escape(self.description)
        esc_default_line_a = cgi.escape(default_line_a)
        esc_default_line_b = cgi.escape(default_line_b)
        lines = [
                esc_description + ':',
                '<br/>',
                _get_textbox_line(esc_label_a, esc_default_line_a, width),
                '<code> -- </code>',
                _get_textbox_line(esc_label_b, esc_default_line_b, width),
                ]
        return lines

    def _validate_interaction(self, first_value, second_value):
        width = second_value - first_value
        if width < 0:
            raise FormError(
                    'the integer in the field "%s" should not be greater than '
                    'the integer in the field "%s"' % (
                        self.first_label, self.second_label))
        if self.low_width is not None:
            if width < self.low_width:
                raise FormError(
                        'the difference between the endpoints '
                        'of the integer interval from '
                        '"%s" to "%s" should be at least %s' % (
                            self.first_label, self.second_label,
                            self.low_width))
        if self.high_width is not None:
            if width > self.high_width:
                raise FormError(
                        'the difference between the endpoints '
                        'of the integer interval from '
                        '"%s" to "%s" should be at most %s' % (
                            self.first_label, self.second_label,
                            self.high_width))

    def process_cmdline_dict(self, d_in, d_out):
        a, b = self._get_temp_items()
        a.process_cmdline_dict(d_in, d_out)
        b.process_cmdline_dict(d_in, d_out)
        self._validate_interaction(d_out[a.label], d_out[b.label])

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object to be decorated with extra attributes
        """
        a, b = self._get_temp_items()
        a.process_fieldstorage(fs)
        b.process_fieldstorage(fs)
        self._validate_interaction(getattr(fs, a.label), getattr(fs, b.label))



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

    def get_preset_pair(self):
        return self.label, self.default_list

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
            raise FormError(
                    'the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError(
                    'the value for the field "%s" is ambiguous' % self.label)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, tuple(self._gen_reduced(value.splitlines())))


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

