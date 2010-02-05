import cgi
from StringIO import StringIO

import numpy

import Util
import MatrixUtil


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
        checked_list = [radio_item for radio_item in radio_items if radio_item.default]
        if len(checked_list) != 1:
            raise FormError('exactly one radio button should be checked by default')
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
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError('no radio button option was selected for the field "%s"' % self.label)
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError('the value for the field "%s" is ambiguous' % self.label)
        # assert that the selected value is actually one of the radio button options
        if value not in set(item.label for item in self.radio_items):
            raise FormError('an invalid radio button option was selected: %s' % value)
        # for the group, set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)
        # for each item, set the value for the attribute in the fieldstorage object
        for radio_item in self.radio_items:
            # verify that the attribute is not already taken
            if hasattr(fs, radio_item.label):
                raise FormError('the object already has the attribute "%s"' % self.label)
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
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # decorate the FieldStorage object with the boolean values of each item separately
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
        @param group_label: the label of the radio button group of which this item is a member
        @return: the list of lines of html text
        """
        lines = []
        # get the escaped label and description
        escaped_label = cgi.escape(self.label)
        escaped_description = cgi.escape(self.description)
        escaped_group_label = cgi.escape(group_label)
        # add the line that makes the checkbox
        radio_button_line_base = 'input type="radio" name="%s" id="%s" value="%s"' % (escaped_group_label, escaped_label, escaped_label)
        if self.default:
            radio_button_line = '<%s checked="yes"/>' % radio_button_line_base
        else:
            radio_button_line = '<%s/>' % radio_button_line_base
        lines.append(radio_button_line)
        # add the line that makes the label
        label_line = '<label for="%s">%s</label>' % (escaped_label, escaped_description)
        lines.append(label_line)
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
        checkbox_line_base = 'input type="checkbox" name="%s" id="%s" value="%s"' % (escaped_label, escaped_label, escaped_label)
        if self.default:
            checkbox_line = '<%s checked="yes"/>' % checkbox_line_base
        else:
            checkbox_line = '<%s/>' % checkbox_line_base
        lines.append(checkbox_line)
        # add the line that makes the label
        label_line = '<label for="%s">%s</label>' % (escaped_label, escaped_description)
        lines.append(label_line)
        # return the list of lines
        return lines

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
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
        label_line = '<label for="%s">%s:</label>' % (esc_label, esc_description)
        lines.append(label_line)
        lines.append('<br/>')
        # add the textbox line
        textbox_line = '<input type="text" name="%s" id="%s" value="%s" size="%d"/>' % (esc_label, esc_label, esc_default_line, width)
        lines.append(textbox_line)
        # return the list of lines
        return lines

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError('the value for the field "%s" is ambiguous' % self.label)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Float:
    """
    An floating point number is requested using a single line in the form.
    """

    def __init__(self, label, description, default_float, low_exclusive=None, low_inclusive=None, high_exclusive=None, high_inclusive=None):
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
        label_line = '<label for="%s">%s:</label>' % (esc_label, esc_description)
        lines.append(label_line)
        lines.append('<br/>')
        # add the textbox line
        textbox_line = '<input type="text" name="%s" id="%s" value="%s" size="%d"/>' % (esc_label, esc_label, esc_default_line, width)
        lines.append(textbox_line)
        # return the list of lines
        return lines

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError('no floating point number was given for the field "%s"' % self.label)
        elif len(values) == 1:
            value_string = values[0]
            try:
                value = float(value_string)
            except ValueError, e:
                raise FormError('%s could not be interpreted as a floating point number' % value_string)
        elif len(values) > 2:
            raise FormError('the value for the field "%s" is ambiguous' % self.label)
        # assert that the floating point number is not out of bounds
        identifier = 'the floating point number in the field "%s"' % self.label
        if self.low_exclusive is not None:
            if value <= self.low_exclusive:
                raise FormError('%s must be greater than %f' % (identifier, self.low_exclusive))
        if self.low_inclusive is not None:
            if value < self.low_inclusive:
                raise FormError('%s must be greater than or equal to %f' % (identifier, self.low_inclusive))
        if self.high_exclusive is not None:
            if value >= self.high_exclusive:
                raise FormError('%s must be less than %f' % (identifier, self.high_exclusive))
        if self.high_inclusive is not None:
            if value > self.high_inclusive:
                raise FormError('%s must be less than or equal to %f' % (identifier, self.high_inclusive))
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Integer:
    """
    An integer is requested using a single line in the form.
    """

    def __init__(self, label, description, default_integer, low=None, high=None):
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
        label_line = '<label for="%s">%s:</label>' % (esc_label, esc_description)
        lines.append(label_line)
        lines.append('<br/>')
        # add the textbox line
        textbox_line = '<input type="text" name="%s" id="%s" value="%s" size="%d"/>' % (esc_label, esc_label, esc_default_line, width)
        lines.append(textbox_line)
        # return the list of lines
        return lines

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError('no integer was given for the field "%s"' % self.label)
        elif len(values) == 1:
            value_string = values[0]
            try:
                value = int(value_string)
            except ValueError, e:
                raise FormError('%s could not be interpreted as an integer' % value_string)
        elif len(values) > 2:
            raise FormError('the value for the field "%s" is ambiguous' % self.label)
        # make sure that the value is in bounds
        if self.low is not None:
            if value < self.low:
                raise FormError('the integer in the field "%s" must be at least %d' % (self.label, self.low))
        if self.high is not None:
            if value > self.high:
                raise FormError('the integer in the field "%s" must be at most %d' % (self.label, self.high))
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)


class Matrix:
    """
    This represents a numpy matrix entered by the user.
    The matrix should be entered in a multi-line text box.
    Each non-empty line should be a row of whitespace separated numbers.
    """

    def __init__(self, label, description, default_matrix, matrix_assertion=None):
        """
        @param label: something like a variable name
        @param description: a single line description of the item
        @param default_matrix: the default numpy array
        @param matrix_assertion: a function of a matrix that may raise MatrixUtil.MatrixError
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
        # get escaped values
        esc_label = cgi.escape(self.label)
        esc_description = cgi.escape(self.description)
        esc_default_string = cgi.escape(MatrixUtil.m_to_string(self.default_matrix))
        # add the label line followed by a line break
        label_line = '<label for="%s">%s:</label>' % (esc_label, esc_description)
        lines.append(label_line)
        lines.append('<br/>')
        # add the textarea header
        textarea_header = '<textarea name="%s" id="%s" rows="%d" cols="70" wrap="off">' % (esc_label, esc_label, nrows)
        lines.append(textarea_header)
        # add the multiple lines of default text
        lines.append(esc_default_string)
        # add the textarea footer
        textarea_footer = '</textarea>'
        lines.append(textarea_footer)
        # return the list of lines
        return lines

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            raise FormError('the value for the field "%s" is empty' % self.label)
        elif len(values) == 1:
            try:
                value = numpy.array(MatrixUtil.read_matrix(StringIO(values[0])))
                if self.matrix_assertion:
                    self.matrix_assertion(value)
            except MatrixUtil.MatrixError, e:
                raise FormError(e)
        elif len(values) > 2:
            raise FormError('the value for the field "%s" is ambiguous' % self.label)
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
        label_line = '<label for="%s">%s:</label>' % (esc_label, esc_description)
        lines.append(label_line)
        lines.append('<br/>')
        # add the textarea header
        textarea_header = '<textarea name="%s" id="%s" rows="%d" cols="70" wrap="off">' % (esc_label, esc_label, nrows)
        lines.append(textarea_header)
        # add the multiple lines of default text
        lines.append(esc_default_string)
        # add the textarea footer
        textarea_footer = '</textarea>'
        lines.append(textarea_footer)
        # return the list of lines
        return lines

    def process_fieldstorage(self, fs):
        """
        @param fs: a FieldStorage object that will be decorated with extra attributes
        """
        # verify that the attribute is not already taken
        if hasattr(fs, self.label):
            raise FormError('the object already has the attribute "%s"' % self.label)
        # read the attribute value from the fieldstorage object
        values = fs.getlist(self.label)
        if not values:
            value = ''
        elif len(values) == 1:
            value = values[0]
        elif len(values) > 2:
            raise FormError('the value for the field "%s" is ambiguous' % self.label)
        # set the value for the attribute in the fieldstorage object
        setattr(fs, self.label, value)

