"""
This is how the output format of a snippet is specified.

Some metaprogramming is done for compatibility with Mobyle.
In Mobyle, all script output is to a file or to stdout.
Also stdout is treated as a file.
In my more constrained system, script output is only to stdout.
"""

from lxml import etree

import Form
import FormOut

def _get_filename_metaprogram(fmt, interpolants):
    if len(interpolants) == 1:
        return '"%s" %% %s' % (fmt, interpolants[0])
    elif len(interpolants) > 1:
        return '"%s" %% (%s)' % (fmt, ', '.join(interpolants))
    else:
        return '"%s"' % fmt

class FormOut(object):

    def __init__(self, filename_format_string, filename_interpolants):
        self.filename_format_string = filename_format_string
        self.filename_interpolants = filename_interpolants

    def get_filename_metaprogram(self):
        """
        This is ugly but necessary metaprogramming for Mobyle compatibility.
        @return: a string that is a python expression
        """
        return _get_filename_metaprogram(
            self.filename_format_string, self.filename_interpolants)

    def get_filename(self, fs):
        """
        @param fs: fieldstorage
        """
        rhs = tuple(getattr(fs, x) for x in self.filename_interpolants)
        return self.filename_format_string % rhs

    def get_contenttype(self, fs):
        return 'text/plain'

    def get_response_headers(self, fs):
        response_headers = []
        response_headers.append(('Content-Type', self.get_contenttype()))
        if hasattr(fs, 'contentdisposition'):
            filename = self.get_filename(fs)
            disposition = '%s; filename=%s' % (fs.contentdisposition, filename)
            response_headers.append(('Content-Disposition', disposition))
        return response_headers

    def get_mobyle_class(self):
        return 'Text'

    def get_mobyle_superclass(self):
        return None

    def _add_mob_xml_datatypes(self, parent):
        etree.SubElement(parent, 'class').text = self.get_mobyle_class()
        superclass = self.get_mobyle_superclass()
        if superclass:
            etree.SubElement(parent, 'superclass').text = superclass

    def add_mob_xml(self, parent, short_name):
        """
        Add a mobyle parameter to the xml tree.
        @param parent: parent etree element
        @param short_name: the short name of the snippet module
        """
        desc = short_name + '_out'
        parameter = etree.SubElement(parent, 'parameter', isout='1')
        etree.SubElement(parameter, 'name').text = desc
        prompt = etree.SubElement(parameter, 'prompt', lang='en')
        prompt.text = desc
        mytype = etree.SubElement(parameter, 'type')
        datatype = etree.SubElement(mytype, 'datatype')
        self._add_mob_xml_datatypes(datatype)
        filenames = etree.SubElement(parameter, 'filenames')
        mycode = etree.SubElement(filenames, 'code', proglang='python')
        mycode.text = self.get_filename_metaprogram()


class Image(FormOut):

    def __init__(self, base_format_string, base_interpolants):
        self.filename_format_string = base_format_string + '.%s'
        self.filename_interpolants = base_interpolants + ['imageformat']

    def get_contenttype(self, fs):
        return Form.g_imageformat_to_contenttype[fs.imageformat]

    def get_mobyle_class(self):
        return 'Picture'

    def get_mobyle_superclass(self):
        return 'Binary'


class Report(FormOut):

    def get_mobyle_class(self):
        return 'Report'

    def get_mobyle_superclass(self):
        return 'Text'


class RTable(Report):
    def get_mobyle_class(self):
        return 'RTable'

class EigenstratInd(Report):
    def get_mobyle_class(self):
        return 'EigenstratInd'

class EigenstratPheno(Report):
    def get_mobyle_class(self):
        return 'EigenstratPheno'

class EigenstratGeno(Report):
    def get_mobyle_class(self):
        return 'EigenstratGeno'

class EigenstratSnp(Report):
    def get_mobyle_class(self):
        return 'EigenstratSnp'
