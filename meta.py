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

class MetaError(Exception): pass

g_common_t1 = [
        'sys', 'os',
        'unittest',
        'itertools', 'csv', 'StringIO']

g_common_t2 = [
        'argparse',
        'pygraphviz', 'pycairo',
        'numpy', 'scipy', 'matplotlib']

g_common_t3 = [
        'Util', 'NewickIO', 'SnippetUtil']

def get_import_paragraphs(raw_lines):
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
            if s.startswith('"""'):
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
    return paragraphs
