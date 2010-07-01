"""
Parse XML files.
"""

import unittest


def indent(elem, level=0):
    """
    @param elem: an xml ElementTree node
    @param level: the current depth
    """
    s = '\n' + '  ' * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = s + "  "
        for e in elem:
            indent(e, level+1)
        if not e.tail or not e.tail.strip():
            e.tail = s
    if level and (not elem.tail or not elem.tail.strip()):
        elem.tail = s


class TestXmlUtil(unittest.TestCase):

    def test_foo(self):
        pass


def main():
    pass

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    #parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', default='-', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False, help='run some unit tests')
    options, args = parser.parse_args()
    # run a test or run a demo
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestXmlUtil)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        main()

