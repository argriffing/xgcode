"""
This module is for picking through the python modules and scripts.

Dependency searching will look at the first few paragraphs
of a python file following the file-level docstring.
The imports should be organized into three paragraphs.
The first paragraph should have imports from the python standard library.
Imports in the second paragraph should be from
external modules, such as are available on pypi.
The third paragraph consists of imports internal to the package;
these are modules from the primordial ooze that do not have
an upstream source.
"""

import collections
import unittest
import itertools
import string
import shutil
import os
import re

import pyparsing
from pyparsing import Word, Literal

import iterutils
import Util

class MetaError(Exception): pass

g_custom_extension_modules = set([
        'Config',
        'day'])

g_common_t1 = set([
        'sys', 'os', 'unittest',
        'collections', 'optparse', 'argparse',
        'unittest', 'math', 'random',
        'itertools', 'csv', 'StringIO'])

g_common_t2 = set([
        'lxml', 'sympy',
        'pygraphviz', 'cairo',
        'numpy', 'scipy', 'matplotlib'])

g_common_t3 = set([
        'smallutil', 'matrixio', 'latexutil',
        'ambignt', 'Progress',
        'EqualArcLayout', 'MatrixUtil', 'iterutils', 'Form', 'Util',
        'Newick', 'FelTree', 'NewickIO', 'SnippetUtil'])

g_default_modules = [
        'auto']

def get_tier(names):
    """
    @param names: a set of module names
    @return: None or the tier
    """
    tiers = set()
    for i, sample in enumerate((g_common_t1, g_common_t2, g_common_t3)):
        if names & sample:
            tiers.add(i)
    if len(tiers) > 1:
        msg_a = 'expected modules from the same tier: '
        msg_b = ', '.join(names)
        raise MetaError(msg_a + msg_b)
    if tiers:
        return list(tiers)[0]

def perm_is_strictly_increasing(perm):
    for a, b in iterutils.pairwise(perm):
        if b <= a:
            return False
    return True

def perm_is_compatible(perm, tiers):
    for v, t in zip(perm, tiers):
        if t is not None:
            if v != t:
                return False
    return True

def _get_import_tiers(paragraphs):
    """
    @param paragraphs: each paragraph has module names
    @return: a tier for each paragraph
    """
    if not paragraphs:
        return tuple()
    paragraphs = [set(p) for p in paragraphs]
    tiers = [get_tier(p) for p in paragraphs]
    valid_perms = []
    for perm in itertools.permutations(range(3), len(paragraphs)):
        if not perm_is_strictly_increasing(perm):
            continue
        if not perm_is_compatible(perm, tiers):
            continue
        valid_perms.append(perm)
    if not valid_perms:
        msg = 'found no valid tier permutation'
        raise MetaError(msg)
    if len(valid_perms) > 1:
        msg = 'found multiple valid pier permutations'
        raise MetaError(msg)
    return valid_perms[0]

def _get_import_paragraphs(raw_lines):
    """
    @param raw_lines: raw lines of a python file
    @return: up to three import paragraphs
    """
    paragraphs = []
    p = []
    if not raw_lines:
        return paragraphs
    # process line by line as a state machine
    doc_start = False
    doc_end = False
    for line in raw_lines:
        if not doc_start:
            s = line.lstrip()
            if not s:
                continue
            if s.startswith('"""') or s.startswith('r"""'):
                doc_start = True
            else:
                msg = 'the file should start with a """ quoted docstring'
                raise MetaError(msg)
            if s.count('"""') > 1:
                doc_end = True
        elif not doc_end:
            s = line.rstrip()
            if not s:
                continue
            if s.endswith('"""'):
                doc_end = True
        else:
            s = line.lstrip()
            if not s:
                if p:
                    paragraphs.append(p)
                    p = []
                continue
            elements = s.split()
            if elements[0] == 'import':
                if ',' in s:
                    msg = 'each import should be on a separate line'
                    raise MetaError(msg)
            if elements[0] in ('import', 'from'):
                p.append(elements[1])
            else:
                break
    if len(paragraphs) > 3:
        raise MetaError('expected at most three import paragraphs')
    return paragraphs

def get_tiered_names(raw_lines):
    """
    @param raw_lines: raw lines of a python file
    @return: three ordered module name sets
    """
    paragraphs = _get_import_paragraphs(raw_lines)
    tiers = _get_import_tiers(paragraphs)
    tiered_names = [set(), set(), set()]
    for t, p in zip(tiers, paragraphs):
        tiered_names[t] = set(p)
    return tiered_names

def _get_const_parser():
    parser = (
            Word(string.letters + '_') +
            Literal('=') +
            Literal('const.read') +
            Literal('(') +
            Literal("'") +
            Word(string.digits + string.lowercase)('dep') +
            Literal("'") +
            Literal(')'))
    return parser

def get_const_deps(raw_lines):
    """
    @param raw_lines: raw lines of a python file
    @return: a set of const-data names
    """
    # look for lines like
    # g_foo = const.read('20100101a')
    deps = set()
    parser = _get_const_parser()
    for line in raw_lines:
        try:
            result = parser.parseString(line)
            if result:
                deps.add(result['dep'])
        except pyparsing.ParseException as e:
            pass
    return deps


class Dep(object):

    def __init__(self):
        self.d_deps = {}

    def get_immediate_deps(self, module_name):
        """
        @param module_name: the name of a module
        @return: three sets of module names
        """
        deps = self.d_deps.get(module_name, None)
        if deps is not None:
            return deps
        filename = module_name + '.py'
        with open(filename) as fin:
            try:
                deps = get_tiered_names(fin)
            except MetaError as e:
                msg = 'dependency format error in %s: %s' % (filename, e)
                raise MetaError(msg)
        self.d_deps[module_name] = deps
        return deps

    def get_transitive_deps(self, module_name):
        deps = [set(), set(), set()]
        visited = set()
        unexpanded = set([module_name]) - set(g_custom_extension_modules)
        while unexpanded:
            name = unexpanded.pop()
            visited.add(name)
            next_deps = self.get_immediate_deps(name)
            deps[0].update(next_deps[0])
            deps[1].update(next_deps[1])
            next_local_deps = next_deps[2]
            next_local_deps -= (visited | set(g_custom_extension_modules))
            deps[2].update(next_local_deps)
            unexpanded.update(next_local_deps)
        return deps

def get_module_and_const_deps(module_names):
    # get transitive module dependencies
    depstate = Dep()
    transitive_deps = [set(), set(), set()]
    for name in module_names:
        filename = name + '.py'
        deps = depstate.get_transitive_deps(name)
        for a, b in zip(transitive_deps, deps):
            a.update(b)
    # get const-data dependencies for all transitive module dependencies
    const_deps = set()
    for name in set(module_names) | transitive_deps[2]:
        filename = name + '.py'
        with open(filename) as fin:
            const_deps.update(get_const_deps(fin))
    # return module and const data deps
    return transitive_deps, const_deps

def _copy_custom_modules(selected_module_names, target):
    selected_filenames = [x + '.py' for x in selected_module_names]
    for name in selected_module_names:
        source_filename = name + '.py'
        target_filename = os.path.join(target, name + '.py')
        shutil.copyfile(source_filename, target_filename)

def _copy_default_modules(target):
    for name in g_default_modules:
        source_filename = name + '.py'
        target_filename = os.path.join(target, name + '.py')
        shutil.copyfile(source_filename, target_filename)

def _copy_const_data(const_deps, target):
    if const_deps:
        const_data_dir = os.path.join(target, 'const-data')
        if not os.path.isdir(const_data_dir):
            os.mkdir(const_data_dir)
        for name in const_deps:
            source_filename = os.path.join('const-data', name + '.dat')
            target_filename = os.path.join(const_data_dir, name + '.dat')
            shutil.copyfile(source_filename, target_filename)

def add_python_files(module_names, python_project):
    module_deps, const_deps = get_module_and_const_deps(module_names)
    selected_module_names = set(module_names) | module_deps[2]
    if not os.path.exists(python_project):
        os.makedirs(python_project)
    _copy_custom_modules(selected_module_names, python_project)
    _copy_default_modules(python_project)
    _copy_const_data(const_deps, python_project)

def get_title(usermod):
    """
    @param usermod: an imported snippet module
    @return: the title of the snippet
    """
    return Util.get_stripped_lines(usermod.__doc__.splitlines())[0]

def get_module_names(manifest, create_all, create_tagged, srcdir='.'):
    """
    This function is for selecting modules.
    @param manifest: None or filename listing modules
    @param create_all: flag
    @param create_tagged: None or a tag
    @return: module names
    """
    pattern = r'^\d{8}[a-zA-Z]\.py$'
    if sum(bool(x) for x in (manifest, create_all, create_tagged)) != 1:
        msg = 'expected exactly one of {manifest, create_all, create_tagged}'
        raise ValueError(msg)
    module_names = []
    if manifest:
        with open(manifest) as fin:
            module_names = [x.strip() for x in fin]
    elif create_all:
        for name in os.listdir(srcdir):
            if re.match(pattern, name):
                module_name = name[:-3]
                module_names.append(module_name)
    elif create_tagged:
        for name in os.listdir(srcdir):
            if re.match(pattern, name):
                module_name = name[:-3]
                usermod = None
                try:
                    usermod = __import__(
                            module_name, globals(), locals(), [], -1)
                except ImportError as e:
                    pass
                if usermod:
                    #FIXME do all modules have g_tags yet?
                    if hasattr(usermod, 'g_tags'):
                        if is_tag_prefix(usermod.g_tags, create_tagged):
                            module_names.append(module_name)
    return module_names

class ImportedModuleInfo:
    """
    Info for one of many imported modules.
    """
    def __init__(self, usermod, identifier):
        """
        @param usermod: the imported module
        @param identifier: a short identifier derived from the description
        """
        self.usermod = usermod
        self.identifier = identifier
    def get_usermod(self):
        return self.usermod
    def get_identifier(self):
        return self.identifier
    def get_name(self):
        return self.usermod.__name__
    def get_title(self):
        return Util.get_stripped_lines(self.usermod.__doc__.splitlines())[0]

def get_usermod_info(module_names, short_name_length):
    """
    @param module_names: generally uninformative names of modules
    @param short_name_length: max length of unique short module names
    @return: a list of module info, and a list of import errors
    """
    import_errors = []
    # The modules need to be imported to get the unique short names.
    usermods = []
    for name in module_names:
        try:
            usermod = __import__(name, globals(), locals(), [], -1)
        except ImportError as e:
            import_errors.append(e)
        usermods.append(usermod)
    # Get all long titles.
    titles = [get_title(usermod) for usermod in usermods]
    # Get corresponding unique short names.
    identifiers = get_short_titles(titles, short_name_length)
    # Get info per module.
    mod_infos = []
    for usermod, identifier in zip(usermods, identifiers):
        mod_infos.append(ImportedModuleInfo(usermod, identifier))
    return mod_infos, import_errors

def is_tag_prefix(tags, prefix):
    """
    The prefix and each tag may be colon-separated.
    @param tags: a list of tags for a module
    @param prefix: look for this prefix tag
    """
    prefix_as_list = prefix.split(':')
    for tag in tags:
        tag_as_list = tag.split(':')
        if Util.list_starts_with(tag_as_list, prefix_as_list):
            return True
    return False

def get_short_titles(titles, length):
    """
    @param titles: a sequence of titles
    @param length: target length of shortened titles
    @return: a sequence of unique shortened titles
    """
    whitelisted = [_get_transformed_title(x, length) for x in titles]
    # if each string occurs only once then we are done
    if len(whitelisted) == len(set(whitelisted)):
        return whitelisted
    # shorten the strings
    maxtrim = len(str(len(titles)))
    whitelisted = [s[:length-maxtrim] for s in whitelisted]
    # get the number of occurrences of each string
    dtotal = collections.defaultdict(int)
    for s in whitelisted:
        dtotal[s] += 1
    # create the new strings
    d = collections.defaultdict(int)
    differentiated = []
    for s in whitelisted:
        if dtotal[s] > 1:
            ntrim = len(str(dtotal[s]))
            s_new = s + str(d[s]+1).zfill(ntrim)
        else:
            s_new = s
        differentiated.append(s_new)
        d[s] += 1
    return differentiated

def _get_transformed_title(title, length):
    """
    @param title: a long and poorly formatted title
    @param length: target length of the shortened title
    @return: a short and less interestingly formatted title
    """
    whitelist = set(string.letters + string.digits)
    # change non-whitelist characters to space
    t = ''.join(c if c in whitelist else ' ' for c in title)
    # remove spaces from the ends and collapse duplicate spaces
    t = '_'.join(t.split())
    # trim to the target length
    t = t[:length]
    # return the lower cased string
    return t.lower()


class TestConstParser(unittest.TestCase):

    def test_const_parser_read(self):
        parser = _get_const_parser()
        line = "foo = const.read('20100101a')"
        result = parser.parseString(line)
        self.assertEqual(result['dep'], '20100101a')

    def test_const_parser_read_rstrip(self):
        parser = _get_const_parser()
        line = "foo = const.read('20100101a').rstrip()"
        result = parser.parseString(line)
        self.assertEqual(result['dep'], '20100101a')


class TestGetShortTitles(unittest.TestCase):

    def test_get_short_titles_a(self):
        titles = [
                '  O   HAI!',
                'this is too long and it is alpha',
                'this is too long and it is beta']
        target_length = 10
        expected = [
                'o_hai',
                'this_is_t1',
                'this_is_t2']
        observed = get_short_titles(titles, target_length)
        self.assertEqual(observed, expected)

    def test_get_short_titles_b(self):
        titles = ['foo']*20
        expected = ['f%02d' % (x+1) for x in range(20)]
        target_length = 3
        observed = get_short_titles(titles, target_length)
        self.assertEqual(observed, expected)

    def test_get_short_titles_tricky(self):
        titles = [
                'asdf1',
                'asdf1',
                'asdf1',
                'asdf2']
        expected = [
                'asdf1',
                'asdf2',
                'asdf3',
                'asdf4']
        target_length = 5
        observed = get_short_titles(titles, target_length)
        self.assertEqual(observed, expected)


if __name__ == '__main__':
    unittest.main()
