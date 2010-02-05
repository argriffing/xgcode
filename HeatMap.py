
import unittest
import sys
from StringIO import StringIO

import Util
import Discretizer
import Monospace


# This sample row major matrix is used to demonstrate matrix visualization.
sample_row_major_matrix = [
        [1, 2, 2, 5, 5, 1, 2, 1, 2, 1],
        [1, 1, 2, 5, 5, 5, 2, 2, 2, 1],
        [1, 2, 2, 2, 5, 6, 5, 1, 2, 1],
        [1, 1, 2, 2, 1, 5, 5, 1, 2, 1],
        [1, 2, 2, 2, 1, 1, 2, 1, 2, 1],
        [1, 2, 2, 2, 2, 1, 2, 2, 2, 1]
        ]

def blue_red_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    assert 0 <= t <= 1
    if t < .5:
        b = 1
        r = 2*t
    else:
        r = 1
        b = 2*(1 - t)
    return (int(r * 255), 0, int(b * 255))

def red_white_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    assert 0 <= t <= 1
    return (255, int(t * 255), int(t * 255))

def blue_white_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    assert 0 <= t <= 1
    return (int(t * 255), int(t * 255), 255)

def black_white_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    assert 0 <= t <= 1
    return (int(t * 255), int(t * 255), int(t * 255))

def white_blue_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    return blue_white_gradient(1 - t)

def white_black_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    return black_white_gradient(1 - t)

def white_red_gradient(t):
    """
    @param t: a fraction between 0 and 1
    @return: an (r, g, b) triple of integers between 0 and 255
    """
    return red_white_gradient(1 - t)



class Legend:

    def __init__(self, values, max_categories, css_prefix, gradient):
        """
        @param values: a set of numeric values
        @param max_categories: the maximum number of categories
        @param css_prefix: the prefix of custom css classes used for background color
        @param gradient: a callable specifying the color gradient
        """
        self.discretizer = Discretizer.Discretizer(values, max_categories)
        self.css_prefix = css_prefix
        self.gradient = gradient

    def index_to_css_class(self, index):
        if index not in range(self.discretizer.get_category_count()):
            raise ValueError('the index is out of range')
        return '%s%d' % (self.css_prefix, index + 1)

    def value_to_css_class(self, value):
        index = self.discretizer.get_category_index(value)
        return self.index_to_css_class(index)

    def gen_style_lines(self):
        ngroups = self.discretizer.get_category_count()
        for i in range(ngroups):
            css_class = self.index_to_css_class(i)
            if ngroups == 1:
                t = 1
            else:
                t = i / float(ngroups - 1)
            yield '.%s {background-color: rgb%s}' % (css_class, self.gradient(t))

    def gen_legend_lines(self):
        for i, (low, high) in enumerate(self.discretizer.get_boundaries()):
            css_class = self.index_to_css_class(i)
            # Compare strings instead of numbers so we do not end up showing
            # [0.333333333333 .. 0.333333333333] for example.
            low_string = '%f' % low
            high_string = '%f' % high
            if low_string == high_string:
                line = '<span class="%s">&nbsp;</span> %s' % (css_class, low_string)
            else:
                line = '<span class="%s">&nbsp;</span> [%s, %s]' % (css_class, low_string, high_string)
            yield line


class HtmlHeatMap:

    def __init__(self, heatmap):
        """
        @param heatmap: a HeatMap object
        """
        self.heatmap = heatmap

    def get_legend(self):
        out = StringIO()
        print >> out, '<code>'
        print >> out, '<br/>\n'.join(self.heatmap.legend.gen_legend_lines())
        print >> out, '</code>'
        return out.getvalue()

    def get_example_html(self):
        """
        @return: the string representing an entire html file
        """
        sio = StringIO()
        print >> sio, '<html>'
        print >> sio, '<head>'
        print >> sio, '<style type="text/css">'
        print >> sio, self.get_style()
        print >> sio, '</style>'
        print >> sio, '</head>'
        print >> sio, '<body>'
        print >> sio, self.get_body()
        print >> sio, '<br/>'
        print >> sio, self.get_legend()
        print >> sio, '</body>'
        print >> sio, '</html>'
        return sio.getvalue()


class TableHeatMap(HtmlHeatMap):

    def get_style(self):
        arr = []
        arr.extend(list(self.heatmap.legend.gen_style_lines()))
        arr.extend(list(self.heatmap._gen_table_style_lines()))
        return '\n'.join(arr)

    def get_body(self):
        arr = []
        arr.append('<table class="heatmap">')
        for row in self.heatmap._gen_table_rows():
            arr.append('<tr>' + row + '</tr>')
        arr.append('</table>')
        return '\n'.join(arr)


class PreHeatMap(HtmlHeatMap):

    def get_style(self):
        return '\n'.join(self.heatmap.legend.gen_style_lines())

    def get_body(self):
        """
        @return: an html pre block with the heatmap data
        """
        arr = []
        arr.append('<pre>')
        arr.extend(list(self.heatmap._gen_pre_lines()))
        arr.append('</pre>')
        return '\n'.join(arr)


class HeatMap:
    """
    Creates a basic heatmap from a row major matrix and a number of heat levels.
    The purpose is simply to visualize a matrix.
    """

    def __init__(self, row_major_matrix, max_categories):
        """
        @param row_major_matrix: a row major matrix of heat values
        @param max_categories: defines the granularity of the discretization
        """
        self.row_major_matrix = row_major_matrix
        values = Util.flattened_nonrecursive(row_major_matrix)
        self.legend = Legend(values, max_categories, 'c', white_black_gradient)

    def _gen_table_style_lines(self):
        """
        These lines go inside the style block in the html head block.
        """
        yield '.heatmap {border: none; border-collapse: collapse; border-spacing: 0}'
        yield '.heatmap td {padding: 0; margin: 0; font-family: monospace;}'

    def _gen_table_rows(self):
        """
        Note that for flexibility the tr tag is not used here.
        @return: an html string defining a row of color blocks within an html table
        """
        for row in self.row_major_matrix:
            arr = []
            for value in row:
                css_class = self.legend.value_to_css_class(value)
                arr.append('<td class="%s">&nbsp;</td>' % css_class)
            yield ''.join(arr)

    def _gen_pre_lines(self):
        """
        Note that for flexibility the pre tag is not used here.
        @return: an html string defining a row of color blocks within an html pre block
        """
        for row in self.row_major_matrix:
            arr = []
            css_classes = [self.legend.value_to_css_class(v) for v in row]
            for css_class, count in Util.rle(css_classes):
                arr.append('<span class="%s">' % css_class)
                arr.append('&nbsp;' * count)
                arr.append('</span>')
            yield ''.join(arr)



class LabeledHeatMap(HeatMap):

    def __init__(self, row_major_matrix, max_categories, row_labels, column_labels):
        """
        @param row_major_matrix: a row major matrix of heat values
        @param max_categories: the maximum number of heat levels
        @param row_labels: short strings labeling the rows
        @param column_labels: short strings labeling the columns
        """
        HeatMap.__init__(self, row_major_matrix, max_categories)
        self.row_labels = row_labels
        self.column_labels = column_labels

    def _get_padded_row_labels(self):
        """
        @return: row labels padded to the length of the longest row label
        """
        rmax = max(len(s) for s in self.row_labels)
        row_labels = [Monospace.left_justify(label, rmax, '.') for label in self.row_labels]
        return row_labels

    def _get_padded_column_labels(self):
        """
        Return a list of padded column labels.
        The number of columns returned is the number of original column labels
        plus the maximum row label length.
        The column labels are padded to equal lengths.
        @return: column labels padded vertically and horizontally
        """
        rmax = max(len(s) for s in self.row_labels)
        cmax = max(len(s) for s in self.column_labels)
        column_labels = [''] * rmax + self.column_labels
        column_labels = [Monospace.left_justify(label, cmax, '.') for label in column_labels]
        return column_labels

    def _gen_table_rows(self):
        """
        For flexibility do not add tr tags.
        """
        row_labels = self._get_padded_row_labels()
        column_labels = self._get_padded_column_labels()
        for row in zip(*column_labels):
            yield ''.join('<td>%s</td>' % c for c in row)
        for label, row_string in zip(row_labels, HeatMap._gen_table_rows(self)):
            yield ''.join('<td>%s</td>' % c for c in label) + row_string

    def _gen_pre_lines(self):
        row_labels = self._get_padded_row_labels()
        column_labels = self._get_padded_column_labels()
        for row in zip(*column_labels):
            yield ''.join(row)
        for label, row_string in zip(row_labels, HeatMap._gen_pre_lines(self)):
            yield label + row_string


class TestHeatMap(unittest.TestCase):

    def test_PreHeatMap(self):
        max_categories = 3
        heatmap = HeatMap(sample_row_major_matrix, max_categories)
        renderer = PreHeatMap(heatmap)
        renderer.get_example_html()

    def test_TableHeatMap(self):
        max_categories = 3
        heatmap = HeatMap(sample_row_major_matrix, max_categories)
        renderer = TableHeatMap(heatmap)
        renderer.get_example_html()

    def test_LabeledHeatMap_pre(self):
        max_categories = 3
        row_labels = list('12345')
        column_labels = list('ABCDEFGHIJ')
        heatmap = LabeledHeatMap(sample_row_major_matrix, max_categories, row_labels, column_labels)
        renderer = PreHeatMap(heatmap)
        renderer.get_example_html()

    def test_LabeledHeatMap_table(self):
        max_categories = 3
        row_labels = list('12345')
        column_labels = list('ABCDEFGHIJ')
        heatmap = LabeledHeatMap(sample_row_major_matrix, max_categories, row_labels, column_labels)
        renderer = TableHeatMap(heatmap)
        renderer.get_example_html()


def run(out, format):
    max_categories = 3
    heatmap = HeatMap(sample_row_major_matrix, max_categories)
    if format == 'pre':
        renderer = PreHeatMap(heatmap)
    elif format == 'table':
        renderer = TableHeatMap(heatmap)
    else:
        raise NotImplementedError('the format "%s" is not implemented' % format)
    out.write(renderer.get_example_html())

def main():
    valid_formats = ('pre', 'table', 'numeric')
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--format', dest='format', default='pre', help='one of ' + str(valid_formats))
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # validate the format option
    if options.format not in valid_formats:
        print 'invalid --format parameter:', options.format, 'is not in', valid_formats
        return
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestHeatMap)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        if options.output_filename == '-':
            out = sys.stdout
        else:
            out = open(options.output_filename, 'w')
        run(out, options.format)
        if out is not sys.stdout:
            out.close()

if __name__ == '__main__':
    main()

