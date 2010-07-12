"""
Mobyle-related code.

The modules Form and FormOut may also contain Mobyle-specific code.
http://infohost.nmt.edu/tcc/help/pubs/pylxml/creating.html
"""

from StringIO import StringIO

from lxml import etree

import Util


class MobyleError(Exception): pass

class FormOutError(MobyleError): pass

def _add_head(auto_path, module_name, doc_lines, parent):
    """
    @param auto_path: path to auto.py
    @param module_name: something like 20100707a
    @param doc_lines: documentation lines
    @param parent: the parent xml level
    """
    head = etree.SubElement(parent, 'head')
    name = etree.SubElement(head, 'name')
    name.text = 'z' + module_name
    version = etree.SubElement(head, 'version')
    version.text = '0.0.1'
    category = etree.SubElement(head, 'category')
    category.text = 'ztools'
    command = etree.SubElement(head, 'command')
    command.text = 'python %s %s' % (auto_path, module_name)
    # add the doc subtree
    doc = etree.SubElement(head, 'doc')
    title = etree.SubElement(doc, 'title')
    title.text = name.text
    description = etree.SubElement(doc, 'description')
    text = etree.SubElement(description, 'text', lang='en')
    text.text = doc_lines[0]

def get_xml(auto_path, module_name):
    """
    @param auto_path: path to auto.py
    @param module_name: something like 20100707a
    """
    # get module info
    usermod = __import__(module_name, globals(), locals(), [], -1)
    form_objects = usermod.get_form()
    if hasattr(usermod, 'get_form_out'):
        form_out = usermod.get_form_out()
    else:
        msg = 'snippet %s provides no output format information' % module_name
        raise FormOutError(msg)
    doc_lines = Util.get_stripped_lines(usermod.__doc__.splitlines())
    # create the root of the tree and the xml document
    program = etree.Element('program')
    doc = etree.ElementTree(program)
    _add_head(auto_path, module_name, doc_lines, program)
    parameters = etree.SubElement(program, 'parameters')
    next_argpos = 1
    for obj in form_objects:
        next_argpos += obj.add_mob_xml(parameters, next_argpos)
    form_out.add_mob_xml(parameters, module_name)
    # serialize the xml
    return etree.tostring(
            doc,
            xml_declaration=True,
            encoding="ISO-8859-1",
            pretty_print=True)

