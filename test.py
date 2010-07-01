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

def gen_module_names():
    """
    Yield names of modules that have been sniffed out.
    """
    for filename in os.listdir('.'):
        if not filename.endswith('.py'):
            continue
        prefix = filename[:-3]
        if not prefix:
            continue
        if not prefix[0].isupper():
            continue
        yield prefix


class TestFinder:
    
    def __init__(self, verbose=False):
        # initialize the important variables
        self.test_classes = []
        self.imported_module_names = []
        self.module_import_errors = []
        # get the list of module names to try
        module_names = list(gen_module_names())
        # create the progress bar
        pbar = Progress.Bar(len(module_names))
        # try to load the modules and get the test classes
        for i, module_name in enumerate(module_names):
            try:
                module = __import__(module_name, globals(), locals())
            except ImportError, e:
                self.module_import_errors.append(e)
            else:
                self.imported_module_names.append(module_name)
                for object_name, object in module.__dict__.items():
                    try:
                        if issubclass(object, unittest.TestCase):
                            self.test_classes.append(object)
                    except TypeError, e:
                        pass
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
    print len(finder.imported_module_names), 'modules were imported successfully'
    print len(finder.module_import_errors), 'modules had import errors'
    # show the import errors
    if finder.module_import_errors:
        print 'import errors:'
        for error in finder.module_import_errors:
            print error

if __name__ == '__main__':
    main()
