"""
Test python snippets.
This script is supposed to test the python snippets
using the default html parameters.
"""

import os
import re
import sys

import numpy as np
import lxml.html as ht

import Form
import Progress

class SnippetTestError(Exception): pass

np.seterr(all='raise')

# failed attempt to catch numpy/scipy warnings
"""
except Warning as e:
    msg = 'warning getting the response: ' + str(e)
    raise SnippetTestError(msg)
"""

# another failed attempt to catch numpy/scipy warnings
"""
with warnings.catch_warnings():
    warnings.simplefilter('error')
"""

def process_module_name(module_name):
    """
    @param module_name: name of the snippet module to try to run
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
    except:
        error_message = str(sys.exc_info()[1])
        msg = 'error importing the snippet: ' + error_message
        raise SnippetTestError(msg)
    try:
        response = module.get_form()
        form_body = Form.get_html_string(response)
        form_html = '<form>' + form_body  + '</form>'
    except:
        error_message = str(sys.exc_info()[1])
        msg = 'get form: ' + error_message
        raise SnippetTestError(msg)
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
        msg = 'default error: ' + str(e)
        raise SnippetTestError(msg)
    # get the result of calling the function
    # using the default parameter values
    if hasattr(module, 'get_response_content'):
        try:
            response = module.get_response_content(mock_field_storage)
            if response is None:
                raise SnippetTestError('no response')
        except:
            error_message = str(sys.exc_info()[1])
            msg = 'get response content: ' + error_message
            raise SnippetTestError(msg)
        else:
            return module, True
    elif hasattr(module, 'get_response'):
        try:
            module.get_response(mock_field_storage)
        except:
            error_message = str(sys.exc_info()[1])
            msg = 'get response: ' + error_message
            raise SnippetTestError(msg)
        else:
            return module, True
    return module, False

def main():
    # first load the module names
    snippet_module_names = []
    for filename in sorted(os.listdir('.')):
        if re.match(r'^\d{8}[a-zA-Z]\.py$', filename):
            prefix = filename.split('.')[0]
            snippet_module_names.append(prefix)
    # report the number of module names detected
    print len(snippet_module_names), 'snippets detected'
    # Try to test each module
    # to assert that no error occurs when the default cgi parameters are used.
    success_count = 0
    for module_name in snippet_module_names:
        print module_name
        module = None
        success = False
        try:
            module, success = process_module_name(module_name)
        except SnippetTestError as e:
            print module_name, ':', e
        if success:
            success_count += 1
        # we are done with the imported snippet
        # so remove it from memory if it was loaded
        if module is not None:
            del module
    print success_count, 'snippets passed without an error'



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
    main()

