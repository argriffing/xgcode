#!/usr/bin/env python

import StringIO
import unittest

import cairo


class CairoUtilError(Exception):
    pass


def hex_string_to_cairo_rgb(hex_string):
    if len(hex_string) != 6:
        raise CairoUtilError('invalid hex color string: %s' % hex_string)
    hex_alphabet = 'abcdefABCDEF0123456789'
    if not (set(hex_string) <= set(hex_alphabet)):
        raise CairoUtilError('invalid hex color string: %s' % hex_string)
    r = int(hex_string[0:2], 16) / 255.0
    g = int(hex_string[2:4], 16) / 255.0
    b = int(hex_string[4:6], 16) / 255.0
    return (r, g, b)

def is_valid_image_format(image_format):
    return image_format in ('png', 'svg', 'ps', 'pdf')

class CairoHelper:
    """
    This class adds image file format flexibility.
    There are two costs for this flexibility:
        - only operations common to all file types are allowed
        - the output image is not streamed but is returned as a string
    """
    def __init__(self, image_format):
        self.image_format = image_format
        self.fout = None
        self.surface = None
        self._validate_image_format()

    def _validate_image_format(self):
        if not self.image_format:
            raise CairoUtilError('no image format was selected')
        if not is_valid_image_format(self.image_format):
            raise CairoUtilError('invalid image format: %s' % self.image_format)

    def create_surface(self, width, height):
        self._validate_image_format()
        if self.fout or self.surface:
            raise CairoUtilError('tried to reuse a surface')
        self.fout = StringIO.StringIO()
        if self.image_format in ('svg', 'png'):
            self.surface = cairo.SVGSurface(self.fout, width, height)
        elif self.image_format == 'ps':
            self.surface = cairo.PSSurface(self.fout, width, height)
        elif self.image_format == 'pdf':
            self.surface = cairo.PDFSurface(self.fout, width, height)
        else:
            raise CairoUtilError('internal image format error')
        return self.surface

    def get_image_string(self):
        """
        After calling this function the surface is invalid.
        """
        self._validate_image_format()
        if not (self.fout or self.surface):
            raise CairoUtilError('no surface is available')
        if self.image_format == 'png':
            out_png = StringIO.StringIO()
            self.surface.write_to_png(out_png)
            self.surface.finish()
            image_string = out_png.getvalue()
        else:
            self.surface.finish()
            image_string = self.fout.getvalue()
        self.image_format = None
        self.fout = None
        self.surface = None
        return image_string


class TestCairoUtil(unittest.TestCase):
    def test_foo(self):
        pass
    def test_bar(self):
        pass


def main():
    """
    Create an example image in each of several formats.
    The example image is from U{cairographics.org/pycairo/}.
    """
    for image_format in ('png', 'svg', 'ps', 'pdf'):
        # set up the drawing surface
        helper = CairoHelper(image_format)
        surface = helper.create_surface(400, 400)
        context = cairo.Context(surface)
        # draw the example image
        context.set_line_width(15)
        context.move_to(200, 100)
        context.line_to(300, 300)
        context.rel_line_to(-200, 0)
        context.close_path()
        context.stroke()
        # write the file
        filename = 'test.' + image_format
        fout = open(filename, 'wb')
        fout.write(helper.get_image_string())
        fout.close()
        print 'created', filename

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestCairoUtil)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()


