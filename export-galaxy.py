"""
Put a project onto a galaxy server.
"""

import re
import os
import sys
import shutil
import subprocess
import tarfile

import argparse
from lxml import etree

import meta
import mobyle
import mobenv
import Util

#FIXME put mobyle functions called here into a non-mobyle and non-galaxy module

def get_name_description_pair(description, name_length=20):
    """
    Galaxy wants the description to be split up.
    The name should be the first part of the input description.
    The description should be the last part of the input description.
    @param description: the long description of the mini app
    @name_length: choose roughly this many letters of the prefix as a name
    """
    elements = description.split()
    k = 0
    for i, x in enumerate(elements):
        if k + len(x) > name_length:
            break
        k += len(x)
    return ' '.join(elements[:i]), ' '.join(elements[i:])

def is_tag_prefix(tags, prefix):
    """
    The prefix and each tag may be colon-separated.
    @param tags: a list of tags for a module
    @param prefix: look for this prefix tag
    """
    #FIXME move this to a non-galaxy-specific module
    prefix_as_list = prefix.split(':')
    for tag in tags:
        tag_as_list = tag.split(':')
        if Util.list_starts_with(tag_as_list, prefix_as_list):
            return True
    return False

def get_module_names(manifest, create_all, create_tagged):
    """
    @param manifest: None or filename listing modules
    @param create_all: flag
    @param create_tagged: None or a tag
    @return: module names
    """
    #FIXME move this to a non-galaxy-specific module
    if sum(bool(x) for x in (manifest, create_all, create_tagged)) != 1:
        msg = 'expected exactly one of {manifest, create_all, create_tagged}'
        raise ValueError(msg)
    module_names = []
    if manifest:
        with open(manifest) as fin:
            module_names = [x.strip() for x in fin]
    elif create_all:
        for name in os.listdir('.'):
            if re.match(r'^\d{8}[a-zA-Z]\.py$', name):
                module_name = name[:-3]
                module_names.append(module_name)
    elif create_tagged:
        for name in os.listdir('.'):
            if re.match(r'^\d{8}[a-zA-Z]\.py$', name):
                module_name = name[:-3]
                usermod = None
                try:
                    usermod = __import__(
                            module_name, globals(), locals(), [], -1)
                except ImportError as e:
                    pass
                if usermod:
                    if hasattr(usermod, 'g_tags'):
                        if is_tag_prefix(usermod.g_tags, create_tagged):
                            module_names.append(module_name)
    return module_names

def make_command(module_name, form_objects):
    """
    The python call is implicit.
    @param module_name: something like '20100707a'
    @param form_objects: a list of input parameter objects
    @return: a single line string with placeholders
    """
    elements = [
            'auto.py', module_name,
            '--output_filename_for_galaxy=$out_file1']
    for obj in form_objects:
        if not obj.web_only():
            elements.append(obj.get_galaxy_cmd())
    return ' '.join(elements)

def get_xml(usermod, module_name, short_name):
    """
    @param usermod: module object
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
    # create the python command string with wildcards
    cmd = make_command(module_name, form_objects)
    # build the xml
    desc_prefix, desc_suffix = get_name_description_pair(doc_lines[0])
    tool = etree.Element('tool', id=short_name, name=desc_prefix)
    if desc_suffix:
        etree.SubElement(tool, 'description').text = desc_suffix
    etree.SubElement(tool, 'command', interpreter='python').text = cmd
    inputs = etree.SubElement(tool, 'inputs')
    outputs = etree.SubElement(tool, 'outputs')
    # add inputs
    for obj in form_objects:
        if not obj.web_only():
            obj.add_galaxy_xml(inputs)
    # add output
    etree.SubElement(outputs, 'data',
            format=form_out.get_galaxy_format(),
            name='out_file1')
    # Add the format tweak if there is an image.
    # This is a hack required because galaxy does not
    # play well with apps which have varying output formats.
    # See the EMBOSS apps for more examples.
    if 'imageformat' in [x.label for x in form_objects]:
        etree.SubElement(tool, 'code', file='galaxy_format_tweak.py')
    # serialize the xml
    return etree.tostring(etree.ElementTree(tool), pretty_print=True)

def add_xml_files(galaxy_root, module_names, short_name_length, tools_subdir):
    """
    @param galaxy_root: root of the galaxy installation
    @param module_names: generally uninformative names of modules
    @param short_name_length: max length of unique short module names
    @param tools_subdir: subdirectory of galaxy_root
    @return: a list of xml filenames, and a list of import errors
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
    short_names = mobyle.get_short_titles(titles, short_name_length)
    # Get the xmls.
    xml_filenames = []
    nsuccesses = 0
    nfailures = 0
    for usermod, name, short_name in zip(usermods, module_names, short_names):
        xml_content = None
        try:
            xml_content = get_xml(usermod, name, short_name)
            nsuccesses += 1
        except:
            error_message = str(sys.exc_info()[1])
            msg = '%s: error making xml: %s' % (name, error_message)
            print >> sys.stderr, msg
            nfailures += 1
        if xml_content:
            xml_filename = short_name + '.xml'
            xml_pathname = os.path.join(
                    galaxy_root, 'tools', tools_subdir, xml_filename)
            with open(xml_pathname, 'w') as fout:
                fout.write(xml_content)
            xml_filenames.append(xml_filename)
    print >> sys.stderr, len(import_errors), 'import errors'
    print >> sys.stderr, nfailures, 'failures to create an xml'
    print >> sys.stderr, nsuccesses, 'successfully created xmls'
    return xml_filenames, import_errors

def get_toolbox_xml(section_name, section_id, tools_subdir, xml_filenames):
    """
    @param section_name: the section name
    @param section_id: the section id
    @param tools_subdir: the name of the tools subdirectory, not a full path
    @param xml_filenames: xml filenames, not full paths
    @return: contents of a toolbox xml file
    """
    # build the xml
    toolbox = etree.Element('toolbox')
    section = etree.SubElement(toolbox, 'section',
            name=section_name, id=section_id)
    # add the items
    for filename in xml_filenames:
        target = os.path.join(tools_subdir, filename)
        etree.SubElement(section, 'tool', file=target)
    # serialize the xml
    return etree.tostring(etree.ElementTree(toolbox), pretty_print=True)

def main(args):
    # get the module names
    module_names = get_module_names(
            args.manifest, args.create_all, args.create_tagged)
    # create the python subtree
    tools_subdir_path = os.path.join(
            args.galaxy_root, 'tools', args.tools_subdir)
    meta.add_python_files(module_names, tools_subdir_path)
    shutil.copyfile('galaxy_format_tweak.py',
            os.path.join(tools_subdir_path, 'galaxy_format_tweak.py'))
    # create the galaxy xml interface files
    xml_filenames, import_errors = add_xml_files(args.galaxy_root,
            module_names, args.short_length, args.tools_subdir)
    for e in import_errors:
        print e
    # create the toolbox xml pointing to the installed xmls
    toolbox_pathname = os.path.join(args.galaxy_root, args.tool_conf)
    section_name = args.tools_subdir
    section_id = args.tools_subdir
    toolbox_xml = get_toolbox_xml(section_name, section_id,
            args.tools_subdir, xml_filenames)
    with open(toolbox_pathname, 'wt') as fout:
        fout.write(toolbox_xml)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',
            help='create xmls for snippets listed in this file')
    parser.add_argument('--create_all', action='store_true',
            help='create xmls for all snippets')
    parser.add_argument('--create_tagged',
            help='create xmls for snippets with this tag')
    parser.add_argument('--galaxy_root', required=True,
            help='path to the Galaxy root directory')
    parser.add_argument('--tool_conf', default='tool_conf.xml.local',
            help='a toolbox xml will be created with this filename')
    parser.add_argument('--tools_subdir', required=True,
            help='python files, const data are copied to this galaxy subdir')
    parser.add_argument('--short_length', type=int, default=20,
            help='max length of shortened snippet names')
    main(parser.parse_args())
