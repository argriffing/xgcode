"""
This script interfaces directly with mod_python and dispatches requests to other scripts.
"""

g_script_directory = '/var/www/python_scripts'
g_extension_directory = '/var/www-extensions'

from StringIO import StringIO
import sys
import os
import re
import cgi

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

def get_doc_root(req):
    doc_root = req.document_root()
    if not doc_root.startswith('/'):
        raise DirectoryError('expected the document root to start with a forward slash: ' + doc_root)
    if not doc_root.endswith('/'):
        # this condition was inconsistent between installs for me
        doc_root += '/'
    return doc_root

def get_script_directory(req):
    """
    @param req: a request object provided by apache2 and mod_python
    @return: the absolute path to the script directory on the local machine
    """
    uri_base = os.path.dirname(req.uri)
    if not uri_base.startswith('/'):
        raise DirectoryError('expected the dirname of req.uri to start with a forward slash')
    doc_root = get_doc_root(req)
    script_directory = doc_root + uri_base[1:]
    whitelist = set()
    whitelist.update(chr(ord('a') + i) for i in range(26))
    whitelist.update(chr(ord('A') + i) for i in range(26))
    whitelist.update(chr(ord('0') + i) for i in range(10))
    whitelist.update(list('_-/'))
    for c in script_directory:
        if c not in whitelist:
            raise DirectoryError('the character "%s" should not exist in the path to the scripts' % c)
    return script_directory

def handler(req):
    """
    If no arguments were provided then display the directory of scripts.
    Otherwise if a valid myscript identifier was provided then dispatch the request to the script.
    Otherwise show an error message.
    """
    import mod_python
    # redirect to the code page if no args were passed
    if not req.args:
        req.content_type = "text/html"
        print >> req, '<head>'
        print >> req, '<meta HTTP-EQUIV="REFRESH" content="0; url=/code">'
        print >> req, '</head>'
        print >> req, '<body></body>'
        return mod_python.apache.OK
    try:
        field_storage = mod_python.util.FieldStorage(req)
        scriptid = field_storage.getfirst('myscript')
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
            # Get the form data.
            form_data = module.get_form()
            # Determine whether we are asking for a form or for a response.
            if len(field_storage.keys()) == 1:
                # the form from the module is either a string or a list of form objects
                if type(form_data) is str:
                    form_html = form_data
                else:
                    form_html = Form.get_html_string(form_data)
                req.content_type = str('text/html')
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
                # possibly parse the field storage data according to the form data
                if type(form_data) is not str:
                    for form_item in form_data:
                        form_item.process_fieldstorage(field_storage)
                # get the response
                content_info, content_text = module.get_response(field_storage)
                for key, value in content_info:
                    if key == 'Content-Type':
                        req.content_type = value
                    else:
                        req.headers_out[key] = value
                req.write(content_text)
        else:
            raise DispatchError('no web interface was found for this script')
    except ImportError as e
        req.content_type = "text/plain"
        print >> req, 'Uncaught ImportError:', e
    except DispatchError as e
        req.content_type = "text/plain"
        print >> req, 'Error:', e
    except HandlingError as e
        req.content_type = "text/plain"
        print >> req, 'Error:', e
    except Form.FormError as e
        req.content_type = "text/plain"
        print >> req, 'Form validation error:', e
    # pretend everything is OK
    return mod_python.apache.OK

def main():
    print 'this is supposed to be called from apache web server'

if __name__ == '__main__':
    main()



