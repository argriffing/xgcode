"""
This module draws a SpatialTree using the cairo API.
This module should only be accessed through its draw() function.
"""

import math
from StringIO import StringIO
from optparse import OptionParser

import cairo

import CairoUtil


def get_tree_image(tree, max_size, image_format):
    """
    Get the image of the tree.
    @param tree: something like a SpatialTree
    @param max_size: (max_width, max_height)
    @param image_format: a string that determines the image format
    @return: a string containing the image data
    """
    # rotate and center the tree on (0, 0)
    tree.fit(max_size)
    # get the width and height of the tree image
    xmin, ymin, xmax, ymax = tree.get_extents()
    width = xmax - xmin
    height = ymax - ymin
    # create the surface
    cairo_helper = CairoUtil.CairoHelper(image_format)
    surface = cairo_helper.create_surface(width, height)
    context = cairo.Context(surface)
    # draw the background
    context.save()
    context.set_source_rgb(.9, .9, .9)
    context.paint()
    context.restore()
    # center on the tree
    context.translate(width/2.0, height/2.0)
    # draw the branches
    for branch in tree.gen_branches():
        context.save()
        color = getattr(branch, 'branch_color', None)
        if color:
            context.set_source_rgb(*CairoUtil.hex_string_to_cairo_rgb(color))
        context.move_to(*branch.src_location)
        context.line_to(*branch.dst_location)
        context.stroke()
        context.restore()
    # get the image string
    return cairo_helper.get_image_string()


class DrawTreeImage:
    """
    This class not currently used.
    This is a class instead of just a namespace so that drawing parameters can be set more flexibly.
    """
    
    def __init__(self, tree, surface):
        self.max_size = (640, 480)
        self.surface = surface
        self.tree = tree

    def draw(self, out):
        """
        @param out: something like a file open for writing
        If the branches are augmented with branch_color members then they will be colored.
        """
        out.write(self._get_rendered_string(self.format, self.max_size, self.tree))

    def _draw_stroke(self, cairo_context, tree):
        """
        Draws the tree defined by the edges in one stroke and finishes the stroke.
        #Draws the tree defined by the edges in one stroke but does not finish the stroke.
        #Calling cairo_context.stroke_extents() after this function returns will get the bounding box of the stroke.
        #Calling cairo_context.stroke() will finish the stroke.
        @param cairo_context: a drawing context
        @param tree: ...
        """
        for branch in tree.gen_branches():
            color = getattr(branch, 'branch_color', None)
            if color:
                cairo_context.set_source_rgb(*CairoUtil.hex_string_to_cairo_rgb(color))
            cairo_context.move_to(*branch.src_location)
            cairo_context.line_to(*branch.dst_location)
            cairo_context.stroke()
            if color:
                cairo_context.set_source_rgb(0, 0, 0)

    def _get_best_sf_theta_extents(self, surface_factory, tree, max_size):
        """
        This function is currently inactive but may be resurrected for more advanced layout.
        @param surface_factory: a function that makes a surface of the required type
        @param tree: ...
        @param max_size: the maximum width and height of the image that would be reasonable for the user to view
        @return: (sf, theta, extents) where theta gives the largest sf scaling factor
        """
        sf_theta_extents_triples = []
        increment_count = 60
        for i in range(increment_count):
            theta = i * (2*math.pi) / increment_count
            """
            # get the extents by doing a dry run
            dummy_file = StringIO()
            dummy_surface = surface_factory(dummy_file, 1, 1)
            dummy_context = cairo.Context(dummy_surface)
            _draw_stroke(dummy_context, list(_gen_rotated_edges(edges, theta)))
            extents = dummy_context.stroke_extents()
            xmin, ymin, xmax, ymax = extents
            dummy_context.stroke()
            dummy_surface.finish()
            dummy_file.close()
            """
            # get the extents analytically
            extents = _compute_extents(list(_gen_rotated_edges(edges, theta)))
            xmin, ymin, xmax, ymax = extents
            # get the width and height of the data
            width = (xmax - xmin)
            height = (ymax - ymin)
            # get the scaling factor that will make the data fit in the box
            current_size = (width, height)
            sf = _get_scaling_factor(current_size, max_size)
            # add to the candidate list
            sf_theta_extents_triples.append((sf, theta, extents))
        return max(sf_theta_extents_triples)

    def _get_rendered_string(self, format, max_size, tree, surface):
        """
        @param format: a string specifying the image format
        @param max_size: the maximum width and height of the image that would be reasonable for the user to view
        @param tree: something like a SpatialTree
        """
        """
        # get normalized edges where the shortest edge has length 10.
        min_edge_length = min(math.hypot(x2-x1, y2-y1) for (x1, y1), (x2, y2) in edges)
        normalized_edges = list(_gen_scaled_edges(edges, 100.0/min_edge_length))
        """
        """
        # get the best parameters by doing dry runs
        sf, theta, (xmin, ymin, xmax, ymax) = _get_best_sf_theta_extents(surface_factory, tree, max_size)
        # get the width and height of the rotated data
        width = (xmax - xmin)
        height = (ymax - ymin)
        current_size = (width, height)
        # define the width and height of the image the user will see
        new_width = sf*width
        new_height = sf*height
        # the points themselves need to be scaled
        # if we use a matrix operation to do the scaling, then line thicknesses and stuff get scaled
        scaled_and_rotated_edges = list(_gen_scaled_edges(_gen_rotated_edges(edges, theta), sf))
        # find the center of the scaled and rotated image
        x_center = sf*float(xmax + xmin) / 2
        y_center = sf*float(ymax + ymin) / 2
        """
        # transform the spatial tree to fit in the box
        tree.fit(max_size)
        #xmin, ymin, xmax, ymax = tree.get_extents()
        #new_width = (xmax - xmin)
        #new_height = (ymax - ymin)
        # emulate a file open for writing
        out = StringIO()
        # create the surface
        #surface = surface_factory(out, new_width, new_height)
        # create the context
        context = cairo.Context(surface)
        # darken the background a little bit
        context.set_source_rgb(.8, .8, .8)
        context.paint()
        context.set_source_rgb(0, 0, 0)
        # center the screen for drawing
        width, height = max_size
        context.translate(width/2.0, height/2.0)
        # draw the edges
        self._draw_stroke(context, tree)
        #context.stroke()
        """
        # extract the data representing the file
        if format == 'png':
            out_png = StringIO()
            surface.write_to_png(out_png)
            surface.finish()
            return out_png.getvalue()
        else:
            surface.finish()
            return out.getvalue()
        """

def main():
    import EqualArcLayout
    import Newick
    import SpatialTree
    # make a spatial tree
    tree = Newick.parse(Newick.daylight_example_tree, SpatialTree.SpatialTree)
    EqualArcLayout.do_layout(tree)
    max_size = (640, 480)
    image_formats = ('png', 'pdf', 'ps', 'svg')
    for image_format in image_formats:
        # get the image string
        image_string = get_tree_image(tree, max_size, image_format)
        # write the image file
        image_filename = 'test.%s' % image_format
        fout = open(image_filename, 'wb')
        fout.write(image_string)
        fout.close()
        print 'created', image_filename
    # write an html file
    html_filename = 'test.html'
    fout = open(html_filename, 'w')
    print >> fout, '<html><body>'
    for image_format in image_formats:
        image_filename = 'test.%s' % image_format
        print >> fout, '<a href="%s">%s</a><br/>' % (image_filename, image_filename)
    print >> fout, '</body></html>'
    fout.close()
    print 'created', html_filename


if __name__ == '__main__':
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    options, args = parser.parse_args()
    main()
