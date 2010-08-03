"""
Code related to the Mobyle environment.
This module should have few dependencies,
because it is used on a target machine as part of the installation.
"""

import os


class EnvironmentInfo:

    def __init__(self, auto_path, python_path, mob_core):
        """
        @param auto_path: path to auto.py
        @param python_path: path to the python executable
        @param mob_core: mobyle core directory
        """
        self.auto_path = auto_path
        self.python_path = python_path
        self.mob_core = mob_core

    def get_deploy_command(self):
        """
        @return: a list to be passed to subprocess.Popen
        """
        return [os.path.join(self.mob_core, 'Tools', 'mobdeploy'), 'deploy']

    def get_xml_dir(self):
        """
        @return: the directory where the xml files will go
        """
        return os.path.join(self.mob_core, 'Local', 'Programs')

    def get_xml_command(self, module_name):
        """
        Get the command string to invoke the snippet without any args.
        @return: the command to be embedded into a Mobyle XML file
        """
        return ' '.join([self.python_path, self.auto_path, module_name])
