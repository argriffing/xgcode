"""
Code related to the Mobyle environment.
This module should have few dependencies,
because it is used on a target machine as part of the installation.
"""

import os


def create_environment_info(auto_path, python_path, mob_core):
    """
    This is a factory function.
    """
    v = get_mobyle_version(mob_core)
    if v == 98:
        klass = EnvironmentInfo98
    elif v == 96:
        klass = EnvironmentInfo96
    else:
        raise ValueError('undetected mobyle version')
    return klass(auto_path, python_path, mob_core)


def get_mobyle_version(mobyle_home):
    set_98 = set(('Doc', 'Example', 'Local',
        'Schema', 'Services', 'Src', 'Tools'))
    if set(os.listdir(mobyle_home)) == set_98:
        return 98
    else:
        return 96


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


class EnvironmentInfo96(EnvironmentInfo):

    def __init__(self, auto_path, python_path, mob_core):
        EnvironmentInfo.__init__(self, auto_path, python_path, mob_core)

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


class EnvironmentInfo98(EnvironmentInfo):

    def __init__(self, auto_path, python_path, mob_core):
        EnvironmentInfo.__init__(self, auto_path, python_path, mob_core)

    def get_xml_dir(self):
        """
        @return: the directory where the xml files will go
        """
        return os.path.join(self.mob_core, 'Local', 'Services', 'Programs')

    def get_xml_command(self, module_name):
        """
        Get the command string to invoke the snippet without any args.
        @return: the command to be embedded into a Mobyle XML file
        """
        bsub = ' '.join((
            'runbsub', '-n 1', '-W 720', '-N', '-q dean',
            '-o result.out', '-e result.err',
            '__bsubArgsEnd__'))
        return ' '.join([bsub, self.python_path, self.auto_path, module_name])
