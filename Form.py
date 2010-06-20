import cgi
from StringIO import StringIO

import numpy as np

import MatrixUtil

g_imageformat_to_contenttype = {
        'svg' : 'image/svg+xml',
        'png' : 'image/png',
        'pdf' : 'application/pdf',
        'postscript' : 'application/postscript'}

g_imageformat_to_ext = {
        'svg' : 'svg',
        'png' : 'png',
        'pdf' : 'pdf',
        'postscript' : 'ps'}

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
        # assert that exactly one radio_item is checked
        checked_list = [item for item in radio_items if item.default]
        if len(checked_list) != 1:
            msg = 'exactly one radio button should be checked by default'
            raise FormError(msg)
        # initialize the member variables
        self.label = label
        self.description = description
        self.radio_items = radio_items

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

    def get_html_lines(self, group_label):
        """
        @param group_label: the label of the radio button group
        @return: the list of lines of html text
        """
        lines = []
        # get the escaped label and description
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_group_label = cgi.escape(group_label)
        # add the line that makes the radiobox
        lines.append(_get_radio_line(esc_group_label, esc_label, self.default))
        # add the line that makes the label
        lines.append(_get_label_line(esc_label, esc_description))
        # return the list of lines
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

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        # get the escaped label and description
        escaped_label = cgi.escape(self.label)
        escaped_description = cgi.escape(self.description)
        # add the line that makes the checkbox
        lines.append(_get_checkbox_line(escaped_label, self.default))
        # add the line that makes the label
        lines.append(_get_label_line(escaped_label, escaped_description))
        # return the list of lines
        return lines

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
        if values:
            is_checked = True
        else:
            is_checked = False
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

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        # calculate a multiple of ten that will hold the string
        width = ((len(self.default_line) / 10) + 1) * 10
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_line = cgi.escape(self.default_line)
        # add the label line followed by a line break
        lines.append(_get_colon_label_line(esc_label, esc_description))
        lines.append('<br/>')
        # add the textbox line
        lines.append(_get_textbox_line(esc_label, esc_default_line, width))
        # return the list of lines
        return lines

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

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        # get the default line of text
        default_line = str(self.default_float)
        # calculate a multiple of ten that will hold the string
        width = ((len(default_line) / 10) + 1) * 10
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_line = cgi.escape(default_line)
        # add the label line followed by a line break
        lines.append(_get_colon_label_line(esc_label, esc_description))
        lines.append('<br/>')
        # add the textbox line
        lines.append(_get_textbox_line(esc_label, esc_default_line, width))
        # return the list of lines
        return lines

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
        # assert that the floating point number is not out of bounds
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

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        # get the default line of text
        default_line = str(self.default_integer)
        # calculate a multiple of ten that will hold the string
        width = ((len(default_line) / 10) + 1) * 10
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_line = cgi.escape(default_line)
        # add the label line followed by a line break
        lines.append(_get_colon_label_line(esc_label, esc_description))
        lines.append('<br/>')
        # add the textbox line
        lines.append(_get_textbox_line(esc_label, esc_default_line, width))
        # return the list of lines
        return lines

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
        # make sure that the value is in bounds
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

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
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
        # add the label line followed by a line break
        lines.append(_get_colon_label_line(esc_label, esc_description))
        lines.append('<br/>')
        # add the textarea header
        lines.append(_get_textarea_header(esc_label, nrows))
        # add the multiple lines of default text
        lines.append(esc_default_string)
        # add the textarea footer
        textarea_footer = '</textarea>'
        lines.append(textarea_footer)
        # return the list of lines
        return lines

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

    def get_html_lines(self):
        """
        @return: the list of lines of html text
        """
        lines = []
        # get the number of rows to use for the textarea
        sio = StringIO(self.default_string)
        nrows = len(list(sio.readlines())) + 1
        nrows = min(nrows, 12)
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_string = cgi.escape(self.default_string)
        # add the label line followed by a line break
        lines.append(_get_colon_label_line(esc_label, esc_description))
        lines.append('<br/>')
        # add the textarea header
        lines.append(_get_textarea_header(esc_label, nrows))
        # add the multiple lines of default text
        lines.append(esc_default_string)
        # add the textarea footer
        textarea_footer = '</textarea>'
        lines.append(textarea_footer)
        # return the list of lines
        return lines

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
            RadioItem('postscript', 'postscript'),
            RadioItem('svg', 'svg')]
        RadioGroup.__init__(self, label, description, radio_items)
