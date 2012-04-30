"""
Test all of the python library files.

Nose will not work directly,
because it requires modules with tests to be named appropriately,
for example by having 'test' in their name or by being in a directory
called 'test'.
"""

import unittest
import os

import Progress

def get_module_names():
    """
    Yield names of modules that have been sniffed out.
    The first pass filters using a file name heuristic.
    The second pass filters using a file content heuristic.
    """
    # every file that has the right filename format passes the first phase
    names_phase_one = []
    for filename in os.listdir('.'):
        if not filename.endswith('.py'):
            continue
        if filename[0].isdigit():
            continue
        if '-' in filename:
            continue
        name = filename[:-3]
        names_phase_one.append(name)
    # files containing all required substrings pass the second phase
    names_phase_two = []
    required_substrings = ('unittest', 'TestCase')
    for name in names_phase_one:
        with open(name + '.py') as fin:
            whole_string = fin.read()
        if all(s in whole_string for s in required_substrings):
            names_phase_two.append(name)
    # return the collection of module names that pass both phases
    return names_phase_two


class TestFinder:
    
    def __init__(self, verbose=False):
        # initialize the variables
        self.test_classes = []
        self.imported_module_names = []
        self.module_import_errors = []
        self.imports_with_tests = []
        # get the list of module names to try
        self.all_module_names = get_module_names()
        # create the progress bar
        pbar = Progress.Bar(len(self.all_module_names))
        # try to load the modules and get the test classes
        for i, module_name in enumerate(self.all_module_names):
            try:
                module = __import__(module_name, globals(), locals())
            except ImportError as e:
                self.module_import_errors.append(e)
            else:
                self.imported_module_names.append(module_name)
                test_classes = []
                for object_name, object in module.__dict__.items():
                    try:
                        if issubclass(object, unittest.TestCase):
                            test_classes.append(object)
                    except TypeError as e:
                        pass
                if test_classes:
                    self.imports_with_tests.append(module_name)
                self.test_classes.extend(test_classes)
            # update the progress bar
            pbar.update(i+1)

def main():
    finder = TestFinder(verbose=True)
    tests = []
    for test_class in finder.test_classes:
        tests.append(unittest.TestLoader().loadTestsFromTestCase(test_class))
    suite = unittest.TestSuite(tests)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # show a summary of the loaded tests
    n_imported = len(finder.imported_module_names)
    n_with_tests = len(finder.imports_with_tests)
    print len(finder.all_module_names), 'potentially testable modules'
    print '  ', n_imported, 'successfully imported'
    print '    ', n_with_tests, 'of these had tests'
    # show the imports without tests
    if n_imported != n_with_tests:
        print 'imported modules that did not have tests:'
        ms = set(finder.imported_module_names) - set(finder.imports_with_tests)
        for name in ms:
            print name
    # show the import errors
    if finder.module_import_errors:
        print 'import errors:'
        for error in finder.module_import_errors:
            print error

if __name__ == '__main__':
    main()
