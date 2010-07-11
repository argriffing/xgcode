"""
This is how the output format of a snippet is specified.

Some metaprogramming is done for compatibility with Mobyle.
In Mobyle, all script output is to a file or to stdout.
Also stdout is treated as a file.
In my more constrained system, script output is only to stdout.
"""

import Form


class FormOut(object):

    def __init__(self, filename_format_string, filename_interpolants):
        self.filename_format_string = filename_format_string
        self.filename_interpolants = filename_interpolants

    def get_filename_code(self):
        """
        This is ugly but necessary metaprogramming for Mobyle compatibility.
        @return: a string that is a python expression
        """
        if len(self.filename_interpolants) == 1:
            return '"%s" %% %s' % (
                    self.filename_format_string,
                    self.filename_interpolants[0])
        elif len(self.filename_interpolants) > 1:
            return '"%s" %% (%s)' % (
                    self.filename_format_string,
                    ', '.join(self.filename_interpolants))
        else:
            return self.filename_format_string

    def get_filename(self, fs):
        """
        @param fs: a fieldstorage object
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


class Image(FormOut):

    def get_contenttype(self, fs):
        return Form.g_imageformat_to_contenttype[fs.imageformat]


class Report(FormOut):
    pass
