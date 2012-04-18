"""
Extra xml formatting.

There is also pretty printing.
import xml.dom.minidom
xml = xml.dom.minidom.parse(xml_fname)
# or xml.dom.minidom.parseString(xml_string)
pretty_xml_as_string = xml.toprettyxml()
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


if __name__ == '__main__':
    unittest.main()

