"""
Test python snippets.

This script is supposed to test the python snippets
using the default parameters.
As of spring 2012 the snippets that use cairo
are having sporadic crashes and hangs when testing using this framework.
These symptoms have not appeared by running the snippets directly on the
command line or through auto.py or through the usual web interface.
"""

import argparse
import os
import re
import sys
import time

import numpy as np
import lxml.html as ht

import Form
import Progress

class SnippetTestError(Exception): pass

np.seterr(all='raise')

def process_module_name(module_name, bad_modules):
    """
    @param module_name: name of the snippet module to try to run
    @param bad_modules: a collection of disallowed module names
    @return: module, success
    """
    module = None
    try:
        requested_methods = [
                '__doc__',
                'get_form',
                'get_response',
                'get_response_content']
        module = __import__(module_name, globals(), locals(),
                requested_methods)
    except Exception as e:
        raise SnippetTestError('error importing the snippet: ' + str(e))
    if any(hasattr(module, x) for x in bad_modules):
        raise SnippetTestError(
                'skipping because the snippet imports a disallowed module')
    try:
        response = module.get_form()
        form_body = Form.get_html_string(response)
        form_html = '<form>' + form_body  + '</form>'
    except Exception as e:
        raise SnippetTestError('get form: ' + str(e))
    # parse the default html parameters from the html string
    document = ht.fragment_fromstring(form_html)
    html_form = document.forms[0]
    html_inputs = html_form.inputs
    # create an object that looks like a FieldStorage object
    mock_field_storage = MockFieldStorage(html_inputs)
    # parse the field storage data according to the form data
    try:
        for form_item in response:
            form_item.process_fieldstorage(mock_field_storage)
    except Form.FormError as e:
        raise SnippetTestError('default error: ' + str(e))
    # get the result of calling the function
    # using the default parameter values
    if hasattr(module, 'get_response_content'):
        try:
            response = module.get_response_content(mock_field_storage)
            if response is None:
                raise SnippetTestError('no response')
        except Exception as e:
            raise SnippetTestError('get response content: ' + str(e))
        else:
            return module, True
    elif hasattr(module, 'get_response'):
        try:
            module.get_response(mock_field_storage)
        except Exception as e:
            raise SnippetTestError('get response: ' + str(e))
        else:
            return module, True
    return module, False

def main(args):
    bad_modules = args.bad_modules.split()
    # first load the module names
    snippet_module_names = []
    for filename in sorted(os.listdir('.')):
        if re.match(r'^\d{8}[a-zA-Z]\.py$', filename):
            prefix = filename.split('.')[0]
            snippet_module_names.append(prefix)
    # Try to test each module
    # to assert that no error occurs when the default cgi parameters are used.
    names_successful = []
    name_to_nseconds = {}
    name_to_linemax = {}
    names_with_unrestored_cwd = []
    for module_name in snippet_module_names:
        print os.getcwd(), module_name
        with open(module_name + '.py') as fin:
            linemax = max(len(line) for line in fin.readlines())
            name_to_linemax[module_name] = linemax
        module = None
        success = False
        t = time.time()
        curdir = os.getcwd()
        try:
            module, success = process_module_name(module_name, bad_modules)
        except SnippetTestError as e:
            print module_name, ':', e
        finally:
            if os.getcwd() != curdir:
                names_with_unrestored_cwd.append(module_name)
            os.chdir(curdir)
        if success:
            nseconds = time.time() - t
            name_to_nseconds[module_name] = nseconds
            names_successful.append(module_name)
        # we are done with the imported snippet
        # so remove it from memory if it was loaded
        if module is not None:
            del module
    print len(snippet_module_names), 'snippets detected'
    print len(names_successful), 'snippets passed without an error'
    print
    # show the snippets with the longest lines
    pairs = [(t, name) for name, t in name_to_linemax.items()]
    bad_pairs = sorted(pairs, reverse=True)[:5]
    print len(bad_pairs), 'snippets with longest lines:'
    for value, name in bad_pairs:
        print name, value
    print
    # show the slowest snippets
    pairs = [(t, name) for name, t in name_to_nseconds.items()]
    bad_pairs = sorted(pairs, reverse=True)[:5]
    print len(bad_pairs), 'slowest snippets with elapsed seconds:'
    for value, name in bad_pairs:
        print name, value
    print
    # show the snippets with unrestored working directory
    if names_with_unrestored_cwd:
        print 'snippets with unrestored working directory:'
        for name in names_with_unrestored_cwd:
            print name
    else:
        print 'all snippets politely restored the working directory'
    print


class MockFieldStorage:
    """
    This is supposed to be like a FieldStorage class.
    """

    def __init__(self, html_inputs):
        """
        @param html_inputs: an lxml FormElement.inputs object
        """
        self.html_inputs = html_inputs

    def getfirst(self, name):
        """
        @param name: the name of the requested parameter
        @return: the first value associated with the requested parameter
        """
        result = self.html_inputs[name]
        value = result.value
        return value

    def getlist(self, name):
        """
        @param name: the name of the requested parameter
        @return: the list of values associated with the requested parameter
        """
        return [self.getfirst(name)]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bad_modules', type=str, default='cairo',
            help='do not test snippets that import these modules')
    main(parser.parse_args())

