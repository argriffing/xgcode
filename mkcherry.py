"""
Copy files as preparation to run the webserver.
"""

import os
import shutil

import CherryUtil


def remove_old_tree():
    shutil.rmtree(CherryUtil.g_root)

def create_xg_directories():
    new_dirs = (
            CherryUtil.g_script_directory,
            CherryUtil.g_epydoc_directory,
            CherryUtil.g_extension_directory)
    for d in new_dirs:
        try:
            os.makedirs(d)
        except OSError:
            pass
        if not os.path.isdir(CherryUtil.g_script_directory):
            raise IOError(d + ' is not a directory')

def copy_xg_files():
    for filename in os.listdir('.'):
        if filename.endswith('.py'):
            pathname = os.path.abspath(filename)
            targetname = os.path.join(CherryUtil.g_script_directory, filename)
            shutil.copyfile(pathname, targetname)

if __name__ == '__main__':
    remove_old_tree()
    create_xg_directories()
    copy_xg_files()
