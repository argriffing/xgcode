"""
Launch a cherrypy server.
"""

from StringIO import StringIO
import sys
import platform
import os
import subprocess
import time
import shutil
import re
import cgi

import argparse
import cherrypy

import SnippetUtil
import Form
import FormHeaderJs

g_web_doc = 'doc'
g_live_doc = 'doc'
g_live_code = 'code'
g_live_log = 'log'

g_script_path = os.path.abspath(sys.argv[0])
g_script_directory = os.path.dirname(g_script_path)

g_current_directory = os.path.abspath(os.curdir)


def get_first_line(multi_line_string):
    """
    @return: a single line string or None
    """
    lines = [line.strip() for line in multi_line_string.splitlines()]
    lines = [line for line in lines if line]
    if lines:
        return lines[0]

class GadgetError(Exception): pass

class Gadget:

    def __init__(self, gadget_path):
        """
        @param gadget_path: the path to the python file of the gadget
        """
        self.module_name = None
        self.module = None
        self.source_link = None
        self.import_error = None
        filename = os.path.basename(gadget_path)
        if not re.match(r'^\d{8}[a-zA-Z]\.py$', filename):
            raise GadgetError('invalid gadget name: ' + filename)
        self.module_name, ext = os.path.splitext(filename)
        try:
            self.module = __import__(self.module_name)
        except ImportError as e:
            self.import_error = e
        epydoc_filename = 'script-' + self.module_name + '_py-pysrc.html'
        if os.path.isfile(os.path.join(g_live_doc, epydoc_filename)):
            self.source_link = '/'.join([g_web_doc, epydoc_filename])

    def __cmp__(self, other):
        return cmp(self.module_name, other.module_name)

    def __str__(self):
        return str([self.module_name, self.module, self.source_link])


def gen_gadgets():
    """
    Yield initialized gadgets.
    Ideally each gadget name is the name of a corresponding module,
    but realistically the named gadgets may be defective in various ways.
    For example, a gadget may fail to import.
    Or a gadget may not have any documentation.
    The gadget directory parameter should be the live code directory.
    @param gadget_directory: the gadget directory
    """
    gadget_directory = g_live_code
    for filename in os.listdir(gadget_directory):
        gadget_path = os.path.join(gadget_directory, filename)
        try:
            yield Gadget(gadget_path)
        except GadgetError as e:
            pass


class FieldStorage(object):

    def __init__(self, param_dict):
        self._param_dict = param_dict

    def getlist(self, key):
        if key in self._param_dict:
            return [self._param_dict[key]]
        else:
            return []

def get_preset_button_tag(preset_index, preset_description):
    return '<button onclick="on_wsf_preset%d();">%s</button>' % (
            preset_index,
            cgi.escape(preset_description))

class GadgetForm(object):

    def __init__(self, module):
        self.module = module
        self.form_objects = None
        self.form_out = None
        self.form_presets = None

    def _init_form(self):
        if self.form_objects is None:
            self.form_objects = self.module.get_form()
        if self.form_out is None:
            self.form_out = self.module.get_form_out()
        if self.form_presets is None:
            if hasattr(self.module, 'get_presets'):
                self.form_presets = self.module.get_presets()
            else:
                self.form_presets = ()

    @cherrypy.expose
    def process(self, **param_dict):
        print param_dict
        start_time = time.time()
        self._init_form()
        fs = FieldStorage(param_dict)
        # determine the content disposition
        if param_dict['submit'] not in ('view', 'download'):
            raise ValueError('invalid content disposition request')
        setattr(fs, 'submit', param_dict['submit'])
        # add the form-specific items to the fs object
        for form_item in self.form_objects:
            form_item.process_fieldstorage(fs)
        # get the header pairs using the new method
        header_pairs = self.form_out.get_response_headers(fs)
        # Use the new get_response_content function if available.
        # Otherwise use the old get_response function.
        if hasattr(self.module, 'get_response_content'):
            content = self.module.get_response_content(fs)
        else:
            deprecated_header_pairs, content = self.module.get_response(fs)
        cherrypy.response.headers.update(dict(header_pairs))
        print 'processed in %f seconds' % (time.time() - start_time)
        return content

    @cherrypy.expose
    def index(self):
        out = StringIO()
        self._init_form()
        form_html = Form.get_html_string(self.form_objects)
        print >> out, '<html>'
        print >> out, '<head>'
        title = SnippetUtil.docstring_to_title(self.module.__doc__)
        if title:
            print >> out, '<title>'
            print >> out, title
            print >> out, '</title>'
        print >> out, FormHeaderJs.get_header_script_text(
                self.form_objects,
                self.form_presets)
        print >> out, '</head>'
        print >> out, '<body>'
        print >> out, SnippetUtil.docstring_to_html(self.module.__doc__)
        print >> out, '<br/><br/>'
        print >> out
        print >> out, '<!-- main form -->'
        print >> out, '<div style="float: left;">'
        print >> out, '<form id="mainform" action="process" method="post">'
        if form_html:
            print >> out, form_html
            print >> out, '<br/><br/>'
        print >> out, '<input type="submit" name="submit" value="view"/>'
        print >> out, '<input type="submit" name="submit" value="download"/>'
        print >> out, '<br/>'
        print >> out, '</form>'
        print >> out, '</div>'
        print >> out
        if self.form_presets:
            print >> out, '<!-- preset configurations -->'
            print >> out, '<div style="float: left;">'
            print >> out, '<fieldset>'
            print >> out, '<legend>preset configurations</legend>'
            button_tags = []
            default_preset = Form.get_default_preset(self.form_objects)
            for i, preset in enumerate([default_preset] + self.form_presets):
                tag = get_preset_button_tag(i, preset.description)
                button_tags.append(tag)
            print >> out, '<br/>\n'.join(button_tags)
            print >> out, '</fieldset>'
            print >> out, '</div>'
            print >> out
        print >> out, '</body>'
        print >> out, '</html>'
        return out.getvalue().rstrip()


class MainForm(object):

    def __init__(self, gadgets):
        self.gadgets = gadgets

    def create_index_contents(self):
        """
        Create the main part of the index html file.
        @return: a chunk of html
        """
        arr = []
        for g in self.gadgets:
            arr.append(g.module_name)
            if g.module:
                link = '[<a href="%s">cgi</a>]' % g.module_name
            else:
                link = '[<span style="color:gray;">cgi</span>]'
            arr.append(link)
            if g.source_link:
                src = '[<a href="%s">src</a>]' % g.source_link
            else:
                src = '[<span style="color:gray;">src</span>]'
            arr.append(src)
            if g.module:
                desc = get_first_line(g.module.__doc__)
                if not desc:
                    desc = '(no description)'
            else:
                desc = '<span style="color:gray;">%s</span>' % g.import_error
            arr.append(desc)
            arr.append('<br />')
        return '\n'.join(arr)

    @cherrypy.expose
    def index(self):
        """
        @return: contents of an html index file.
        """
        arr = []
        start_time = time.time()
        arr += ['<html><body>']
        arr += ['<code>']
        arr += [self.create_index_contents()]
        arr += ['</code>']
        nseconds = time.time() - start_time
        arr += ['<!-- generated in %.2f seconds -->' % nseconds]
        arr += ['</body></html>']
        return '\n'.join(arr)

def get_static_conf():
    current_directory = os.path.abspath(os.curdir)
    doc_directory = os.path.join(current_directory, g_live_doc)
    conf = {'/doc': {
        'tools.staticdir.on': True,
        'tools.staticdir.dir': doc_directory}}
    return conf

def gen_extension_paths():
    extensions_dir = os.path.join(g_script_directory, 'extensions')
    for d in os.listdir(extensions_dir):
        extension_path = os.path.join(extensions_dir, d)
        if os.path.isdir(extension_path):
            yield extension_path

def build_extensions():
    extension_paths = list(gen_extension_paths())
    for ext in extension_paths:
        cmd = ('python', 'setup.py', 'build_ext')
        proc = subprocess.Popen(cmd, cwd=ext, stdout=subprocess.PIPE)
        output = proc.communicate()[0]
        log_filename = os.path.join(g_live_log, os.path.basename(ext) + '.log')
        with open(log_filename, 'wt') as fout:
            fout.write(output)
        build_dir = os.path.join(ext, 'build')
        dlib = '-'.join([
            '.'.join(['lib', platform.uname()[0].lower()]),
            platform.processor(),
            '.'.join(platform.python_version_tuple()[:2])])
        library_dir = os.path.join(build_dir, dlib)
        for sofile in os.listdir(library_dir):
            sofile_path = os.path.join(library_dir, sofile)
            if sofile.endswith('.so') and os.path.isfile(sofile_path):
                sofile_target = os.path.join(g_live_code, sofile)
                shutil.copyfile(sofile_path, sofile_target)

def gen_module_paths(source_dir):
    """
    @param source_dir: a directory
    """
    for f in os.listdir(source_dir):
        fpath = os.path.join(source_dir, f)
        if os.path.isfile(fpath):
            if f.endswith('.so') or f.endswith('.py'):
                yield fpath

def create_documentation():
    module_paths = list(gen_module_paths(g_live_code))
    cmd = ['epydoc', '--output=' + g_live_doc] + module_paths
    output = subprocess.Popen(cmd, stdout=subprocess.PIPE).communicate()[0]
    log_filename = os.path.join(g_live_log, 'epydoc.log')
    with open(log_filename, 'wt') as fout:
        fout.write(output)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--host', type=str, default='127.0.0.1',
            help='use 0.0.0.0 for all public ports')
    parser.add_argument('--port', type=int, default=8080,
            help='should be at least 1024 unless you have root access')
    parser.add_argument('--mkdocs', action='store_true', help='build docs')
    args = parser.parse_args()
    if g_current_directory == g_script_directory:
        raise ValueError('Run this script from a temporary "live" directory.')
    if os.path.isdir(g_live_code):
        shutil.rmtree(g_live_code)
    os.makedirs(g_live_code)
    if os.path.isdir(g_live_log):
        shutil.rmtree(g_live_log)
    os.makedirs(g_live_log)
    for filename in os.listdir(g_script_directory):
        if filename.endswith('.py'):
            src = os.path.abspath(os.path.join(g_script_directory, filename))
            dst = os.path.abspath(os.path.join(g_live_code, filename))
            shutil.copyfile(src, dst)
    const_data_src = os.path.join(g_script_directory, 'const-data')
    const_data_dst = 'code/const-data'
    if os.path.isdir(const_data_dst):
        shutil.rmtree(const_data_dst)
    shutil.copytree(const_data_src, const_data_dst)
    build_extensions()
    if args.mkdocs:
        if os.path.isdir(g_live_doc):
            shutil.rmtree(g_live_doc)
        create_documentation()
    sys.path.remove(g_script_directory)
    sys.path.append(os.path.abspath(g_live_code))
    gadgets = list(reversed(sorted(gen_gadgets())))
    cherrypy.config.update({
        'server.socket_host': args.host,
        'server.socket_port': args.port})
    main_form = MainForm(gadgets)
    for g in gadgets:
        if g.module:
            form = GadgetForm(g.module)
            setattr(main_form, g.module_name, form)
    cherrypy.quickstart(main_form, '/', config=get_static_conf())

