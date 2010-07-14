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

from lxml import etree

import Util
import iterutils


class MobyleError(Exception): pass

class FormOutError(MobyleError): pass

def _get_transformed_title(title, length):
    """
    @param title: a long and poorly formatted title
    @param length: target length of the shortened title
    @return: a short and less interestingly formatted title
    """
    whitelist = set(string.letters + string.digits)
    # change non-whitelist characters to space
    t = ''.join(c if c in whitelist else ' ' for c in title)
    # remove spaces from the ends and collapse duplicate spaces
    t = '_'.join(t.split())
    # trim to the target length
    t = t[:length]
    # return the lower cased string
    return t.lower()

def get_short_titles(titles, length):
    """
    @param titles: a sequence of titles
    @param length: target length of shortened titles
    @return: a sequence of unique shortened titles
    """
    whitelisted = [_get_transformed_title(x, length) for x in titles]
    # if each string occurs only once then we are done
    if len(whitelisted) == len(set(whitelisted)):
        return whitelisted
    # shorten the strings
    maxtrim = len(str(len(titles)))
    whitelisted = [s[:length-maxtrim] for s in whitelisted]
    # get the number of occurrences of each string
    dtotal = collections.defaultdict(int)
    for s in whitelisted:
        dtotal[s] += 1
    # create the new strings
    d = collections.defaultdict(int)
    differentiated = []
    for s in whitelisted:
        if dtotal[s] > 1:
            ntrim = len(str(dtotal[s]))
            s_new = s + str(d[s]+1).zfill(ntrim)
        else:
            s_new = s
        differentiated.append(s_new)
        d[s] += 1
    return differentiated

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

def get_xml(usermod, auto_path, module_name, short_name):
    """
    @param usermod: module object
    @param auto_path: path to auto.py
    @param module_name: something like '20100707a'
    @param short_name: something like 'plot_pca_3d'
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
    etree.SubElement(head, 'category').text = 'all'
    # add categories for input
    input_categories = [obj.__class__.__name__ for obj in form_objects]
    for input_category in iterutils.unique_everseen(input_categories):
        etree.SubElement(head, 'category').text = 'input:' + input_category
    # add categories for output
    output_category = form_out.__class__.__name__
    etree.SubElement(head, 'category').text = 'output:' + output_category
    # add categories for tags
    for tag in tags:
        etree.SubElement(head, 'category').text = tag
    command = etree.SubElement(head, 'command')
    command.text = 'python %s %s' % (auto_path, module_name)
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

def add_xml_files(module_names, path_to_auto, xml_target, short_name_length):
    """
    @param module_names: generally uninformative names of modules
    @param path_to_auto: a path that will go into the xml
    @param xml_target: the xml files will go into this path
    @param short_name_length: max length of unique short module names
    @return: a list of import errors
    """
    import_errors = []
    # The modules need to be imported to get the unique short names.
    usermods = []
    for name in module_names:
        try:
            usermod = __import__(name, globals(), locals(), [], -1)
        except ImportError as e:
            import_errors.append(e)
        usermods.append(usermod)
    # Get all long titles.
    titles = []
    for usermod in usermods:
        doc_lines = Util.get_stripped_lines(usermod.__doc__.splitlines())
        titles.append(doc_lines[0])
    # Get corresponding unique short names.
    short_names = get_short_titles(titles, short_name_length)
    for usermod, name, short_name in zip(usermods, module_names, short_names):
        try:
            xml_content = get_xml(usermod, path_to_auto, name, short_name)
        except MobyleError as e:
            xml_content = None
            print e
        if xml_content:
            xml_filename = os.path.join(xml_target, short_name + '.xml')
            with open(xml_filename, 'w') as fout:
                fout.write(xml_content)
    return import_errors

class TestMobyle(unittest.TestCase):

    def test_get_short_titles_a(self):
        titles = [
                '  O   HAI!',
                'this is too long and it is alpha',
                'this is too long and it is beta']
        target_length = 10
        expected = [
                'o_hai',
                'this_is_t1',
                'this_is_t2']
        observed = get_short_titles(titles, target_length)
        self.assertEqual(observed, expected)

    def test_get_short_titles_b(self):
        titles = ['foo']*20
        expected = ['f%02d' % (x+1) for x in range(20)]
        target_length = 3
        observed = get_short_titles(titles, target_length)
        self.assertEqual(observed, expected)

    def test_get_short_titles_tricky(self):
        titles = [
                'asdf1',
                'asdf1',
                'asdf1',
                'asdf2']
        expected = [
                'asdf1',
                'asdf2',
                'asdf3',
                'asdf4']
        target_length = 5
        observed = get_short_titles(titles, target_length)
        self.assertEqual(observed, expected)

if __name__ == '__main__':
    unittest.main()
