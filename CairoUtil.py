"""
Tools for pycairo.
"""

from StringIO import StringIO
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
        self.fout = StringIO()
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
            out_png = StringIO()
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

def create_test_image(image_format, width, height):
    """
    @param image_format: a string like an image extension
    @param width: width in pixels
    @param height: height in pixels
    @return: the contents of the image as a string
    """
    # do some initialization
    helper = CairoHelper(image_format)
    surface = helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw an off-white background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # get the contents of the image
    return helper.get_image_string()


class TestCairoUtil(unittest.TestCase):

    def test_create_test_image(self):
        width = 640
        height = 480
        for image_format in ('png', 'svg', 'ps', 'pdf'):
            create_test_image(image_format, width, height)


def create_demo_files():
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
    unittest.main()
