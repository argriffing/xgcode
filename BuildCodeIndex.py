"""
Create a static html page
which is an index of the web executable python scripts.
"""

g_script_directory = '.'
g_extension_directory = '.'
g_doc_directory = 'web-build/epydoc-html'

g_web_script_directory = '/python_scripts'
g_web_doc_directory = '/phydoc'

import StringIO
import sys
import os
import re
import cgi
import time

def add_to_path(directory):
    if directory not in sys.path:
        sys.path.append(directory)

add_to_path(g_script_directory)
add_to_path(g_extension_directory)

from SnippetUtil import HandlingError
import SnippetUtil
import Form

class DispatchError(Exception): pass
class DirectoryError(Exception): pass

def write_directory_html(script_directory, doc_directory, out):
    """
    Display a directory of available scripts.
    @param script_directory: the absolute path to the script directory
    @param out: an object that has the same interface as an open output file
    """
    start_time = time.time()
    print >> out, '<html><body>'
    print >> out, '<code>'
    write_directory_html_guts(script_directory, doc_directory, out)
    print >> out, '</code>'
    stop_time = time.time()
    print >> out, '<!--'
    print >> out, 'generated in', stop_time - start_time, 'seconds'
    print >> out, '-->'
    print >> out, '</body></html>'

def write_directory_html_guts(script_directory, doc_directory, out):
    """
    Display a directory of available scripts.
    @param script_directory: the absolute path to the script directory
    @param out: an object that has the same interface as an open output file
    """
    for filename in reversed(sorted(os.listdir(script_directory))):
        if not re.match(r'^\d{8}[a-zA-Z]\.py$', filename):
            continue
        prefix = filename.split('.')[0]
        try:
            doc_filename = 'script-%s_py-pysrc.html' % prefix
            doc_path = os.path.join(doc_directory, doc_filename)
            web_doc_path = os.path.join(g_web_doc_directory, doc_filename)
            web_script_path = os.path.join(g_web_script_directory, 'myscript.py')
            module = __import__(prefix, globals(), locals(), ['__doc__', 'handler', 'get_form', 'get_response'])
            action_link = '<a href="%s?myscript=%s">cgi</a>' % (web_script_path, prefix)
            if os.path.isfile(doc_path):
                doc_link = '<a href="%s">src</a>' % web_doc_path
            else:
                doc_link = '<span style="color:gray;">src</span>'
            module_description_lines = [line for line in module.__doc__.split('\n') if line.strip()]
            if module_description_lines:
                description = module_description_lines[0].strip()
            else:
                description = '(no description)'
            print >> out, '%s [%s] [%s] %s<br/>' % (prefix, action_link, doc_link, description)
            del module
        except ImportError, e:
            action_link = '<span style="color:gray;">cgi</span>'
            doc_link = '<span style="color:gray;">src</span>'
            description = '<span style="color:gray;">%s</span>' % cgi.escape(str(e))
            print >> out, '%s [%s] [%s] %s<br/>' % (prefix, action_link, doc_link, description)

def main():
    """
    Write an html file.
    """
    write_directory_html(g_script_directory, g_doc_directory, sys.stdout)

if __name__ == '__main__':
    main()


