"""
Code related to the Mobyle environment.
This module should have few dependencies,
because it is used on a target machine as part of the installation.
"""

import os


def create_environment_info(auto_path, python_path, mob_core, mob_version):
    """
    This is a factory function.
    @param auto_path: path to auto.py
    @param python_path: path to the python executable
    @param mob_core: mobyle core directory
    @param mob_version: mobyle version as a string
    """
    if mob_version == '0.98':
        klass = EnvironmentInfo98
    elif mob_version == '0.96':
        klass = EnvironmentInfo96
    else:
        raise ValueError('unknown mobyle version ' + mob_version)
    return klass(auto_path, python_path, mob_core)

class EnvironmentInfo:

    def __init__(self, auto_path, python_path, mob_core):
        self.auto_path = auto_path
        self.python_path = python_path
        self.mob_core = mob_core

    def get_index_command(self):
        """
        @return: a list to be passed to subprocess.Popen
        """
        return [os.path.join(self.mob_core, 'Tools', 'mobdeploy'), 'index']

    def get_clean_command(self):
        """
        @return: a list to be passed to subprocess.Popen
        """
        return [os.path.join(self.mob_core, 'Tools', 'mobdeploy'), 'clean']

    def get_deploy_command(self):
        """
        @return: a list to be passed to subprocess.Popen
        """
        return [os.path.join(self.mob_core, 'Tools', 'mobdeploy'), 'deploy']

    def get_xml_command(self, module_name):
        """
        Get the command string to invoke the snippet without any args.
        @return: the command to be embedded into a Mobyle XML file
        """
        return ' '.join([self.python_path, self.auto_path, module_name])


class EnvironmentInfo96(EnvironmentInfo):

    def __init__(self, auto_path, python_path, mob_core):
        EnvironmentInfo.__init__(self, auto_path, python_path, mob_core)

    def get_xml_dir(self):
        """
        @return: the directory where the xml files will go
        """
        return os.path.join(self.mob_core, 'Local', 'Programs')


class EnvironmentInfo98(EnvironmentInfo):

    def __init__(self, auto_path, python_path, mob_core):
        EnvironmentInfo.__init__(self, auto_path, python_path, mob_core)

    def get_xml_dir(self):
        """
        @return: the directory where the xml files will go
        """
        return os.path.join(self.mob_core, 'Local', 'Services', 'Programs')
