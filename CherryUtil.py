import os
import sys

def add_to_path(directory):
    if not os.path.isdir(directory):
        raise ValueError(directory + ' is not a valid directory')
    if directory not in sys.path:
        sys.path.append(directory)

try:
    add_to_path(os.path.expanduser('~/.xgconfig'))
except ValueError, e:
    msg = ' '.join(['The configuration directory was not found:', e])
    raise ValueError(msg)

import xgconfig

g_root = xgconfig.g_root_directory
g_script_directory = os.path.join(g_root, 'script')
g_extension_directory = os.path.join(g_root, 'extension')
g_epydoc_directory = os.path.join(g_root, 'static', 'epydoc')
