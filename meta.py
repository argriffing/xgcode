"""This is the main application for wsgi dispatching.
"""

import StringIO
import sys
import os
import re
import cgi


script_directory = '/var/www/wsgi'
extension_directory = '/var/www-extensions'
doc_directory = '/var/www/phydoc'

web_path_to_scripts = '/w'
web_path_to_docs = '/phydoc'

def add_to_path(directory):
    if directory not in sys.path:
        sys.path.append(directory)

add_to_path(script_directory)
add_to_path(extension_directory)

import SnippetUtil

def scriptid_is_valid(scriptid):
    """
    @param scriptid: should be something like 20080401c
    @return: None when the scriptid does not match the expected format
    """
    pattern = r'\d{8}[a-z]'
    return re.match(pattern, scriptid)

def application(environ, start_response):
    """
    This application function is called by the wsgi framework.
    @param environ: a dictionary of environment variables
    @param start_response: a function taking (status, response_headers) parameters.
    @return: the html page or results data
    """
    out = StringIO.StringIO()
    status = '200 OK'
    response_headers = [('Content-Type', 'text/plain')]
    try:
        fs = cgi.FieldStorage(fp=environ['wsgi.input'], environ=environ)
        scriptid = fs.getfirst('scriptid')
        if scriptid is None:
            # no snippet id was specified so show the listing
            response_headers = [('Content-Type', 'text/html')]
            print >> out, get_directory_html(script_directory, doc_directory)
        elif scriptid_is_valid(scriptid):
            # get a snippet form or a snippet response
            module = __import__(scriptid, globals(), locals(), ['get_form', 'get_response'])
            if 'getform' in fs:
                response_headers = [('Content-Type', 'text/html')]
                form_text = module.get_form()
                print >> out, '<html>'
                print >> out, '<body>'
                print >> out, SnippetUtil.docstring_to_html(module.__doc__)
                print >> out, '<br/><br/>'
                print >> out, '<form method="post">'
                print >> out, '<input type="hidden" name="scriptid" value="%s"/>' % scriptid
                print >> out, '<input type="hidden" name="getresponse" value="1"/>'
                print >> out, form_text
                print >> out, '<br/><br/>'
                print >> out, '<input type="submit" name="mysubmit" value="Submit"/><br/>'
                print >> out, '</form>'
                print >> out, '</body>'
                print >> out, '</html>'
            elif 'getresponse' in fs:
                response_headers, text = module.get_response(fs)
                print >> out, text
        else:
            # there was an error
            for key, value in environ.items():
                print >> out, key, ':', value
            print >> out, '--- FieldStorage ---'
            for key in fs:
                value = fs[key]
                print >> out, key, ':', value
    except Exception, e:
        response_headers = [('Content-Type', 'text/plain')]
        print >> out, 'Fail:'
        print >> out, e
    start_response(status, response_headers)
    return [out.getvalue().strip()]

def get_directory_html(script_directory, doc_directory):
    """
    Display a directory of available scripts.
    @param script_directory: the absolute path to the script directory on the local machine
    @param doc_directory: the path to the documentation
    """
    out = StringIO.StringIO()
    print >> out, '<html><body>'
    print >> out, '<code>'
    for filename in sorted(os.listdir(script_directory)):
        if not re.match(r'^\d{8}[a-zA-Z]\.py$', filename):
            continue
        prefix = filename.split('.')[0]
        try:
            doc_filename = 'script-%s_py-pysrc.html' % prefix
            doc_path = os.path.join(doc_directory, doc_filename)
            expected_function_names = ['get_form', 'get_response']
            module = __import__(prefix, globals(), locals(), ['__doc__'] + expected_function_names)
            # define the cgi link
            cgi_link = '<a href="%s/meta.py?scriptid=%s&getform=1">cgi</a>' % (web_path_to_scripts, prefix)
            for attr in expected_function_names:
                if not hasattr(module, attr):
                    cgi_link = '<span style="color:gray;">cgi</span>'
            # define the src link
            if os.path.isfile(doc_path):
                src_link = '<a href="%s/%s">src</a>' % (web_path_to_docs, doc_filename)
            else:
                src_link = '<span style="color:gray;">src</span>'
            module_description_lines = [line for line in module.__doc__.split('\n') if line.strip()]
            if module_description_lines:
                description = module_description_lines[0].strip()
            else:
                description = '(no description)'
            print >> out, '%s [%s] [%s] %s<br/>' % (prefix, cgi_link, src_link, description)
            del module
        except ImportError:
            print >> out, '%s (Error)' % prefix
    print >> out, '</code>'
    print >> out, '</body></html>'
    return out.getvalue().strip()

def handler(req):
    """
    This handler function is called by the mod_python framework.
    If no arguments were provided then display the directory of scripts.
    Otherwise if a valid myscript identifier was provided then dispatch the request to the script.
    Otherwise show an error message.
    @param req: a mod_python request object
    """
    import mod_python

if __name__ == '__main__':
    print get_directory_html(script_directory, doc_directory)
