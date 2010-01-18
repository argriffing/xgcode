"""
Create a Markov model of English language text.
This is for testing various Markov model algorithms.
"""

import unittest
import os
import StringIO
import logging
import sys

import Util
import Config

# a directory readable and writable from the web
data_directory = Config.data_path

# a module-wide logging object
logger = logging.getLogger('EnglishModel')

# these letters are the possible states in the simple text model
simple_text_states = 'abcdefghijklmnopqrstuvwxyz '

# do not show debugging messages when running as an imported module
logger.setLevel(logging.ERROR)
logger.addHandler(logging.StreamHandler(sys.stdout))


def get_transition_matrix():
    """
    @return: a transition matrix in convenient dictionary form
    """
    count_matrix = _get_count_matrix()
    transition_matrix = {}
    for a in simple_text_states:
        row_sum = float(sum(count_matrix[(a, b)] for b in simple_text_states))
        for b in simple_text_states:
            pair = (a, b)
            transition_matrix[pair] = count_matrix[pair] / row_sum
    return transition_matrix

def _get_count_matrix():
    """
    Get the count matrix with caching.
    @return: a count matrix in convenient dictionary form
    """
    file_path = os.path.join(data_directory, 'tale_of_two_cities_counts.dat')
    if os.path.isfile(file_path):
        logger.debug('found the cached count matrix file')
        fin = open(file_path, 'r')
        count_matrix = _read_count_matrix(fin)
        fin.close()
    else:
        logger.debug('failed to find the cached count matrix file')
        raw_text = _get_raw_text()
        sio = StringIO.StringIO(raw_text)
        simple_text = _raw_text_to_simple_text(sio)
        logger.debug('processed the raw text into simple text')
        count_matrix = dict(((a, b), 1) for a in simple_text_states for b in simple_text_states)
        for a, b in zip(simple_text[:-1], simple_text[1:]):
            count_matrix[(a, b)] += 1
        logger.debug('created the count matrix')
        fout = open(file_path, 'w')
        _write_count_matrix(fout, count_matrix)
        fout.close()
        logger.debug('cached the count matrix')
    return count_matrix

def _write_count_matrix(fout, cm):
    """
    @param fout: an open file for writing
    @param cm: a count matrix in convenient dictionary form
    """
    for (a, b), count in sorted(cm.items()):
        print >> fout, '%s%s:%d' % (a, b, count)

def _read_count_matrix(fin):
    """
    @param fin: an open file for reading
    @return: a count matrix in convenient dictionary form
    """
    cm = {}
    for line in fin:
        if line.count(':') != 1:
            continue
        transition_str, count_str = line.split(':')
        assert len(transition_str) == 2, transition_str
        a, b = transition_str
        cm[(a, b)] = int(count_str)
    return cm


def _raw_text_to_simple_text(fin):
    """
    @param fin: an ascii text file open for reading
    @return: a simplified text string
    """
    arr = []
    space = True
    for c in fin.read():
        if c == "'":
            continue
        elif c.isalpha():
            if space:
                arr.append(' ')
            arr.append(c.lower())
            space = False
        else:
            space = True
    if space:
        arr.append(' ')
    return ''.join(arr)

def _get_raw_text():
    """
    Get some text with caching.
    @return: the text string of Tale of Two Cities.
    """
    file_path = os.path.join(data_directory, 'tale_of_two_cities.txt')
    if os.path.isfile(file_path):
        logger.debug('found the cached raw text file')
        fin = open(file_path, 'r')
        text_string = fin.read()
        fin.close()
        logger.debug('read %d characters from the cached file' % len(text_string))
    else:
        logger.debug('failed to find the cached raw text file')
        text_string = _download_raw_text()
        logger.debug('read %d characters from project gutenberg' % len(text_string))
        fout = open(file_path, 'w')
        fout.write(text_string)
        fout.close()
        logger.debug('cached the raw text file')
    return text_string

def _download_raw_text():
    """
    @return: the text string of Tale of Two Cities.
    """
    import urllib
    url = 'http://www.gutenberg.org/files/98/98.txt'
    furl = urllib.urlopen(url)
    text_string = furl.read()
    furl.close()
    return text_string

class TestEnglishModel(unittest.TestCase):
    def test_foo(self):
        pass
    def test_bar(self):
        pass

def main():
    logger.debug('hello')
    transition_matrix = get_transition_matrix()
    logger.debug('got a transition matrix with %d transitions' % len(transition_matrix))
    logger.debug('goodbye')

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose', default=False)
    #parser.add_option('-o', '--output', dest='output_filename', metavar='FILE', help='output file')
    parser.add_option('--test', action='store_true', dest='test', default=False)
    options, args = parser.parse_args()
    if options.test:
        suite = unittest.TestLoader().loadTestsFromTestCase(TestEnglishModel)
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        if options.verbose:
            logger.setLevel(logging.DEBUG)
        main()


