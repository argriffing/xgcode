"""
This experimental script will run other scripts.
It should import the script as a module and then
use the get_request function in the module to
get the interface.
"""

import sys

import argparse

import Form

if __name__ == '__main__':
    usage = 'example usage:  $ python cmd.py 20100623a'
    if len(sys.argv) != 2:
        raise ValueError(usage)
    script_name, module_name = sys.argv
    usermod = __import__(module_name, globals(), locals(), [], -1)
    form_objects = usermod.get_form()
    print Form.get_help_string(form_objects)
