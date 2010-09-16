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

class ImportedModuleInfo:
    """
    Info for one of many imported modules.
    """
    def __init__(self, usermod, identifier):
        """
        @param usermod: the imported module
        @param identifier: a short identifier derived from the description
        """
        self.usermod = usermod
        self.identifier = identifier
    def get_usermod(self):
        return self.usermod
    def get_identifier(self):
        return self.identifier
    def get_name(self):
        return self.usermod.__name__
    def get_title(self):
        return Util.get_stripped_lines(self.usermod.__doc__.splitlines())[0]

def get_split_title(description, prefix_length=20):
    """
    Galaxy wants the description to be split up.
    The name should be the first part of the input description.
    The description should be the last part of the input description.
    @param description: the long description of the mini app
    @prefix_length: choose roughly this many letters of the prefix as a name
    """
    elements = description.split()
    k = 0
    for i, x in enumerate(elements):
        if k + len(x) > prefix_length:
            break
        k += len(x)
    return ' '.join(elements[:i]), ' '.join(elements[i:])

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
    Get the galaxy XML describing a single tool.
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
    desc_prefix, desc_suffix = get_split_title(doc_lines[0])
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

def get_usermod_info(module_names, short_name_length):
    """
    @param module_names: generally uninformative names of modules
    @param short_name_length: max length of unique short module names
    @return: a list of module info, and a list of import errors
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
    titles = [meta.get_title(usermod) for usermod in usermods]
    # Get corresponding unique short names.
    identifiers = mobyle.get_short_titles(titles, short_name_length)
    # Get info per module.
    mod_infos = []
    for usermod, identifier in zip(usermods, identifiers):
        mod_infos.append(ImportedModuleInfo(usermod, identifier))
    return mod_infos, import_errors

def add_xml_files(galaxy_root, module_names, short_name_length, tools_subdir):
    """
    The XML files are added under two conditions.
    Under the first condition the files are added directly to galaxy.
    In this case only the XML filenames are subsequently needed to
    create the toolbox XML file.
    Under the second condition the files are added to a tool suite archive.
    In this case the id, name, and description of the added modules
    are needed to create the suite configuration XML file.
    This function is concerned with the first case.
    @param galaxy_root: root of the galaxy installation
    @param module_names: generally uninformative names of modules
    @param short_name_length: max length of unique short module names
    @param tools_subdir: subdirectory of galaxy_root
    @return: a list of xml filenames, and a list of import errors
    """
    mod_infos, import_errors = get_usermod_info(module_names, short_name_length)
    # Get the xmls.
    xml_filenames = []
    nsuccesses = 0
    nfailures = 0
    for info in mod_infos:
        usermod = info.get_usermod()
        name = info.get_name()
        short_name = info.get_identifier()
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

def add_xml_archive_files(module_names, short_name_length, archive):
    """
    The XML files are added under two conditions.
    Under the first condition the files are added directly to galaxy.
    In this case only the XML filenames are subsequently needed to
    create the toolbox XML file.
    Under the second condition the files are added to a tool suite archive.
    In this case the id, name, and description of the added modules
    are needed to create the suite configuration XML file.
    This function is concerned with the second case.
    @param module_names: generally uninformative names of modules
    @param short_name_length: max length of unique short module names
    @param archive: the path to the output directory to be tarred
    @return: a list of added module infos, and a list of import errors
    """
    mod_infos, import_errors = get_usermod_info(module_names, short_name_length)
    # Get the xmls.
    added_infos = []
    nsuccesses = 0
    nfailures = 0
    for info in mod_infos:
        usermod = info.get_usermod()
        name = info.get_name()
        short_name = info.get_identifier()
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
            xml_pathname = os.path.join(archive, short_name + '.xml')
            with open(xml_pathname, 'w') as fout:
                fout.write(xml_content)
            added_infos.append(info)
    print >> sys.stderr, len(import_errors), 'import errors'
    print >> sys.stderr, nfailures, 'failures to create an xml'
    print >> sys.stderr, nsuccesses, 'successfully created xmls'
    return added_infos, import_errors

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

def get_suite_config_xml(added_infos, suite_name):
    """
    Create suite_config.xml contents for a galaxy tool suite archive.
    @param added_infos: ImportedModuleInfo objects for added tools
    @param suite_name: the name of the suite directory, not the whole path
    @return: contents of suite_config.xml
    """
    suite = etree.Element('suite', id=suite_name,
            name='Suite of misc tools', version='1.0.0')
    suite_description = 'Suite of misc tools for Galaxy'
    etree.SubElement(suite, 'description').text = suite_description
    for info in added_infos:
        module_id = info.get_identifier()
        module_name = get_split_title(info.get_title())[0]
        module_description = info.get_title()
        tool = etree.SubElement(suite, 'tool',
                id=module_id, name=module_name, version='1.0.0')
        etree.SubElement(tool, 'description').text = module_description
    return etree.tostring(etree.ElementTree(suite), pretty_print=True)

def main_non_archive(args):
    # validation
    if not args.galaxy_root:
        msg = 'in non-archive mode the galaxy root must be specified'
        raise ValueError(msg)
    if not args.tools_subdir:
        msg = 'in non-archive mode the tools subdirectory must be specified'
        raise ValueError(msg)
    # get the module names
    module_names = meta.get_module_names(
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

def main_archive(args):
    # validation
    if args.galaxy_root:
        msg = 'in archive mode the galaxy root must not be specified'
        raise ValueError(msg)
    if args.tools_subdir:
        msg = 'in archive mode the tools subdirectory must not be specified'
        raise ValueError(msg)
    # define the archive extension and the compression command
    archive_extension = '.tar.bz2'
    archive_prefix = os.path.basename(args.suite_archive.rstrip('/'))
    archive_name = archive_prefix + archive_extension
    archive_cmd = ['tar', 'cjvf', archive_name, args.suite_archive]
    # delete the suite directory and archive if they exist
    try:
        shutil.rmtree(args.suite_archive)
        os.remove(archive_name)
    except OSError as e:
        pass
    # create the empty suite directory
    os.makedirs(args.suite_archive)
    # get the module names
    module_names = meta.get_module_names(
            args.manifest, args.create_all, args.create_tagged)
    # add the python files
    meta.add_python_files(module_names, args.suite_archive)
    shutil.copyfile('galaxy_format_tweak.py',
            os.path.join(args.suite_archive, 'galaxy_format_tweak.py'))
    # create the galaxy xml interface files
    mod_infos, import_errors = add_xml_archive_files(
            module_names, args.short_length, args.suite_archive)
    for e in import_errors:
        print e
    # create the toolbox xml pointing to the installed xmls
    config_pathname = os.path.join(args.suite_archive, 'suite_config.xml')
    config_xml = get_suite_config_xml(mod_infos, archive_prefix)
    with open(config_pathname, 'wt') as fout:
        fout.write(config_xml)
    # use subprocess instead of tarfile to create the tgz
    subprocess.call(archive_cmd)

def main(args):
    if args.suite_archive:
        main_archive(args)
    else:
        main_non_archive(args)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',
            help='create xmls for snippets listed in this file')
    parser.add_argument('--create_all', action='store_true',
            help='create xmls for all snippets')
    parser.add_argument('--create_tagged',
            help='create xmls for snippets with this tag')
    parser.add_argument('--short_length', type=int, default=20,
            help='max length of shortened snippet names')
    parser.add_argument('--suite_archive',
            help='create this directory and a corresponding archive')
    parser.add_argument('--galaxy_root',
            help='path to the Galaxy root directory')
    parser.add_argument('--tool_conf', default='tool_conf.xml.local',
            help='a toolbox xml will be created with this filename')
    parser.add_argument('--tools_subdir',
            help='python files, const data are copied to this galaxy subdir')
    main(parser.parse_args())
