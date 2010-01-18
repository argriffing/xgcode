#! /usr/bin/env python

"""
This script acts as a cgi adapter for my python snippets.
The cgi calls are processes instead of threads,
so that it might be possible to use SIGALRM to limit the time spent in a cgi script.
This is not possible in mod_python for example,
because each request is in its own thread instead of its own process.
"""

import StringIO
import sys
import os
import re
import cgi
import cgitb
import signal

cgitb.enable()

script_directory = '/var/www/python_scripts'
extension_directory = '/var/www-extensions'
doc_directory = '/var/www/phydoc'

def add_to_path(directory):
    if directory not in sys.path:
        sys.path.append(directory)

add_to_path(script_directory)
add_to_path(extension_directory)

from SnippetUtil import HandlingError
import SnippetUtil

class DispatchError(Exception): pass
class DirectoryError(Exception): pass

def scriptid_is_valid(scriptid):
    """
    @param scriptid: should be something like 20080401c
    @return: None when the scriptid does not match the expected format
    """
    pattern = r'\d{8}[a-z]'
    return re.match(pattern, scriptid)

def get_doc_root(req):
    doc_root = req.document_root()
    if not doc_root.startswith('/'):
        raise DirectoryError('expected the document root to start with a forward slash')
    if not doc_root.endswith('/'):
        raise DirectoryError('expected the document root to end with a forward slash')
    return doc_root

def get_script_directory(script_filename):
    """
    @param script_filename: the full path to the script filename; also SCRIPT_FILENAME passed by cgi
    @return: the absolute path to the script directory on the local machine
    """
    """
    uri_base = os.path.dirname(req.uri)
    if not uri_base.startswith('/'):
        raise DirectoryError('expected the dirname of req.uri to start with a forward slash')
    doc_root = get_doc_root(req)
    script_directory = doc_root + uri_base[1:]
    """
    #script_directory = os.path.dirname(script_filename)
    whitelist = set()
    whitelist.update(chr(ord('a') + i) for i in range(26))
    whitelist.update(chr(ord('A') + i) for i in range(26))
    whitelist.update(chr(ord('0') + i) for i in range(10))
    whitelist.update(list('_-/'))
    for c in script_directory:
        if c not in whitelist:
            raise DirectoryError('the character "%s" should not exist in the path to the scripts' % c)
    return script_directory

def write_directory_html(script_directory, doc_directory, out):
    """
    Display a directory of available scripts.
    @param script_directory: the absolute path to the script directory on the local machine
    @param out: an object that has the same interface as an open output file
    """
    print >> out, '<html><body>'
    print >> out, '<code>'
    for filename in reversed(sorted(os.listdir(script_directory))):
        if not re.match(r'^\d{8}[a-zA-Z]\.py$', filename):
            continue
        prefix = filename.split('.')[0]
        try:
            doc_filename = 'script-%s_py-pysrc.html' % prefix
            doc_path = os.path.join(doc_directory, doc_filename)
            module = __import__(prefix, globals(), locals(), ['__doc__', 'handler', 'get_form', 'get_response'])
            action_link = '<a href="main.py?myscript=%s">cgi</a>' % prefix
            if os.path.isfile(doc_path):
                doc_link = '<a href="/phydoc/%s">src</a>' % doc_filename
            else:
                doc_link = '<span style="color:gray;">src</span>'
            module_description_lines = [line for line in module.__doc__.split('\n') if line.strip()]
            if module_description_lines:
                description = module_description_lines[0].strip()
            else:
                description = '(no description)'
            print >> out, '%s [%s] [%s] %s<br/>' % (prefix, action_link, doc_link, description)
            del module
        except ImportError:
            print >> out, '%s [import error]<br/>' % prefix
    print >> out, '</code>'
    print >> out, '</body></html>'

def handler(req):
    """
    If no arguments were provided then display the directory of scripts.
    Otherwise if a valid myscript identifier was provided then dispatch the request to the script.
    Otherwise show an error message.
    """
    try:
        field_storage = mod_python.util.FieldStorage(req)
        scriptid = field_storage.getfirst('myscript')
        if not req.args:
            # if the page was called with no arguments then show a directory html page
            try:
                script_directory = get_script_directory(req)
                doc_directory = os.path.join(get_doc_root(req), 'phydoc')
            except DirectoryError, e:
                req.content_type = "text/plain"
                print >> req, 'Error:', e
            else:
                req.content_type = "text/html"
                page_buffer = StringIO.StringIO()
                write_directory_html(script_directory, doc_directory, page_buffer)
                req.write(page_buffer.getvalue())
        else:
            # the page was called with arguments so get a response from the module
            if not scriptid:
                raise DispatchError('no script snippet was specified')
            if not re.match(r'^\d{8}[a-zA-Z]$', scriptid):
                raise DispatchError('the specified script name did not meet the format requirements')
            try:
                module = __import__(scriptid, globals(), locals(), ['__doc__', 'handler', 'get_form', 'get_response'])
            except ImportError:
                raise DispatchError('the script could not be imported')
            # process the module differently depending on its type
            if hasattr(module, 'get_form') and hasattr(module, 'get_response'):
                # This is the more advanced dispatch method.
                # Determine whether we are asking for a form or for a response.
                if len(field_storage.keys()) == 1:
                    # send a form
                    form_html = module.get_form()
                    req.content_type = "text/html"
                    print >> req, '<html>'
                    title = SnippetUtil.docstring_to_title(module.__doc__)
                    if title:
                        print >> req, '<head>'
                        print >> req, '<title>'
                        print >> req, title
                        print >> req, '</title>'
                        print >> req, '</head>'
                    print >> req, '<body>'
                    print >> req, SnippetUtil.docstring_to_html(module.__doc__)
                    print >> req, '<br/><br/>'
                    print >> req, '<form method="post">'
                    print >> req, '<input type="hidden" name="myscript" value="%s"/>' % scriptid
                    print >> req, form_html
                    print >> req, '<br/><br/>'
                    print >> req, '<input type="submit" name="mysubmit" value="Submit"/><br/>'
                    print >> req, '</form>'
                    print >> req, '</body>'
                    print >> req, '</html>'
                else:
                    content_info, content_text = module.get_response(field_storage)
                    for key, value in content_info:
                        if key == 'Content-Type':
                            req.content_type = value
                        else:
                            req.headers_out[key] = value
                    req.write(content_text)
            else:
                raise DispatchError('no web interface was found for this script')
    except DispatchError, e:
        req.content_type = "text/plain"
        print >> req, 'Error:', e
    except HandlingError, e:
        req.content_type = "text/plain"
        print >> req, 'Error:', e
    # pretend everything is OK
    return mod_python.apache.OK

def do_cgi():
    """
    This is called when the script is run as a cgi script.
    """
    # get the fields sent by the browser
    fs = cgi.FieldStorage()
    # start writing a response
    out = StringIO.StringIO()
    # initialize the header dictionary
    header_dict = {}
    # get some environment variables from the cgi
    document_root = os.getenv('DOCUMENT_ROOT')
    script_filename = os.getenv('SCRIPT_FILENAME')
    # try to generate some useful data to send to the user
    try:
        if not list(fs):
            # if the page was called with no arguments then show a directory html page
            doc_directory = os.path.join(document_root, 'phydoc')
            header_dict['Content-Type'] = 'text/html'
            write_directory_html(script_directory, doc_directory, out)
        else:
            # the page was called with arguments so get a response from the module
            scriptid = fs.getfirst('myscript')
            if not scriptid:
                raise DispatchError('no script snippet was specified')
            if not re.match(r'^\d{8}[a-zA-Z]$', scriptid):
                raise DispatchError('the specified script name did not meet the format requirements')
            try:
                module = __import__(scriptid, globals(), locals(), ['__doc__', 'handler', 'get_form', 'get_response'])
            except ImportError:
                raise DispatchError('the script could not be imported')
            # process the module differently depending on its type
            if hasattr(module, 'get_form') and hasattr(module, 'get_response'):
                # This is the more advanced dispatch method.
                # Determine whether we are asking for a form or for a response.
                if len(fs.keys()) == 1:
                    # send a form
                    form_html = module.get_form()
                    header_dict['Content-Type'] = 'text/html'
                    print >> out, '<html>'
                    title = SnippetUtil.docstring_to_title(module.__doc__)
                    if title:
                        print >> out, '<head>'
                        print >> out, '<title>'
                        print >> out, title
                        print >> out, '</title>'
                        print >> out, '</head>'
                    print >> out, '<body>'
                    print >> out, SnippetUtil.docstring_to_html(module.__doc__)
                    print >> out, '<br/><br/>'
                    print >> out, '<form method="post">'
                    print >> out, '<input type="hidden" name="myscript" value="%s"/>' % scriptid
                    print >> out, form_html
                    print >> out, '<br/><br/>'
                    print >> out, '<input type="submit" name="mysubmit" value="Submit"/><br/>'
                    print >> out, '</form>'
                    print >> out, '</body>'
                    print >> out, '</html>'
                else:
                    content_info, content_text = module.get_response(fs)
                    header_dict.update(content_info)
                    out.write(content_text)
            else:
                raise DispatchError('no web interface was found for this script')
    except (DirectoryError, DispatchError, HandlingError), e:
        header_dict['Content-Type'] = 'text/plain'
        print >> out, 'Error:', e
    # write the headers
    for key, value in header_dict.items():
        print key + ':', value
    # write a blank line after the headers
    print
    # write the data
    sys.stdout.write(out.getvalue())

if __name__ == '__main__':
    # assume that the script is called as a cgi script
#    print 'Content-Type: text/plain'
#    print
#    print 'hello world'
#    vars = ('DOCUMENT_NAME', 'DOCUMENT_ROOT', 'DOCUMENT_URI', 'PATH_INFO')
#    for var in vars:
#        print var + ':', os.getenv(var)
#    cgi.print_environ()
    do_cgi()
