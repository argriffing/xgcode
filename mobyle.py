"""
Mobyle-related code.

The modules Form and FormOut may also contain Mobyle-specific code.
http://infohost.nmt.edu/tcc/help/pubs/pylxml/creating.html
"""

from StringIO import StringIO
import os
import unittest
import string
import collections
import sys

from lxml import etree

import Util
import iterutils
import meta


class MobyleError(Exception): pass

class FormOutError(MobyleError): pass

def _add_redirection_parameter(parent, next_argpos):
    """
    Add a mobyle parameter to the xml tree.
    @param parent: parent etree element
    @param next_argpos: a 1-based integer for cmdline arg ordering
    @return: the number of args added on the command line
    """
    mobyle_class = 'Boolean'
    meta_code = '" --write_to_file_for_mobyle"'
    parameter = etree.SubElement(parent, 'parameter', ishidden='1')
    etree.SubElement(parameter, 'name').text = 'write_to_file_for_mobyle'
    prompt = etree.SubElement(parameter, 'prompt', lang='en')
    prompt.text = 'write to file for mobyle'
    mytype = etree.SubElement(parameter, 'type')
    datatype = etree.SubElement(mytype, 'datatype')
    etree.SubElement(datatype, 'class').text = mobyle_class
    fmt = etree.SubElement(parameter, 'format')
    code = etree.SubElement(fmt, 'code', proglang='python')
    code.text = meta_code
    etree.SubElement(parameter, 'argpos').text = '%d' % next_argpos
    return 1


class CategoryInfo:
    """
    Specify how the mobyle xmls should be categorized.
    """
    def __init__(self, show_io_types, show_tags, universal):
        """
        @param show_io_types: True if we should categorize by io type
        @param show_tags: True if we should categorize by tag
        @param universal: None or a category encompassing all xmls
        """
        self.show_io_types = show_io_types
        self.show_tags = show_tags
        self.universal = universal

    def gen_categories(self, tags, form_objects, form_out):
        """
        @return: an iterable of category strings
        """
        if self.show_io_types:
            cats = [obj.__class__.__name__ for obj in form_objects]
            for cat in iterutils.unique_everseen(cats):
                yield 'input:' + cat
            if form_out:
                yield 'output:' + form_out.__class__.__name__
        if self.show_tags:
            for tag in tags:
                yield tag
        if self.universal:
            yield self.universal



def get_xml(cat_info, env_info, usermod, module_name, short_name, runbsub):
    """
    @param cat_info: how xmls will be categorized
    @param env_info: relevant places in the filesystem
    @param usermod: module object
    @param module_name: something like '20100707a'
    @param short_name: something like 'plot_pca_3d'
    @param runbsub: a command
    """
    # get module info
    form_objects = usermod.get_form()
    if hasattr(usermod, 'get_form_out'):
        form_out = usermod.get_form_out()
    else:
        msg = 'snippet %s provides no output format information' % module_name
        raise FormOutError(msg)
    doc_lines = Util.get_stripped_lines(usermod.__doc__.splitlines())
    try:
        tags = usermod.g_tags
    except AttributeError:
        tags = []
    # create the root of the tree and the xml document
    program = etree.Element('program')
    doc = etree.ElementTree(program)
    # add the head subtree
    head = etree.SubElement(program, 'head')
    etree.SubElement(head, 'name').text = short_name
    etree.SubElement(head, 'version').text = '0.0.1'
    # add categories
    for cat in cat_info.gen_categories(tags, form_objects, form_out):
        etree.SubElement(head, 'category').text = cat
    command = etree.SubElement(head, 'command')
    cmd_text = env_info.get_xml_command(module_name)
    if runbsub:
        cmd_text = ' '.join((runbsub, cmd_text))
    command.text = cmd_text
    # add the head.doc subtree
    subtree_doc = etree.SubElement(head, 'doc')
    title = etree.SubElement(subtree_doc, 'title').text = short_name
    description = etree.SubElement(subtree_doc, 'description')
    etree.SubElement(description, 'text', lang='en').text = doc_lines[0]
    # add the parameters
    parameters = etree.SubElement(program, 'parameters')
    next_argpos = 1
    for obj in form_objects:
        if not obj.web_only():
            next_argpos += obj.add_mob_xml(parameters, next_argpos)
    next_argpos += _add_redirection_parameter(parameters, next_argpos)
    form_out.add_mob_xml(parameters, short_name)
    # serialize the xml
    return etree.tostring(
            doc,
            xml_declaration=True,
            encoding="ISO-8859-1",
            pretty_print=True)

def add_xml_files(cat_info, env_info,
        module_names, short_name_length, local_xml_dir, runbsub):
    """
    @param cat_info: how xmls will be categorized
    @param env_info: relevant places in the filesystem
    @param module_names: generally uninformative names of modules
    @param short_name_length: max length of unique short module names
    @param local_xml_dir: usually the xml dir in env_info
    #param runbsub: a path to a python script
    @return: a command
    """
    mod_infos, import_errors = meta.get_usermod_info(
            module_names, short_name_length)
    # Get the xmls.
    nsuccesses = 0
    nfailures = 0
    for info in mod_infos:
        usermod = info.get_usermod()
        name = info.get_name()
        short_name = info.get_identifier()
        xml_content = None
        try:
            xml_content = get_xml(
                    cat_info, env_info, usermod, name, short_name, runbsub)
            nsuccesses += 1
        except:
            error_message = str(sys.exc_info()[1])
            msg = '%s: error making xml: %s' % (name, error_message)
            print >> sys.stderr, msg
            nfailures += 1
        if xml_content:
            xml_filename = os.path.join(local_xml_dir, short_name + '.xml')
            with open(xml_filename, 'w') as fout:
                fout.write(xml_content)
    print >> sys.stderr, len(import_errors), 'import errors'
    print >> sys.stderr, nfailures, 'failures to create an xml'
    print >> sys.stderr, nsuccesses, 'successfully created xmls'
    return [x.get_identifier() for x in mod_infos], import_errors

class TestMobyle(unittest.TestCase):
    pass
if __name__ == '__main__':
    unittest.main()
