"""Do some SNP analysis. [UNFINISHED]
"""

from StringIO import StringIO
import time
import itertools
import optparse
import random

from SnippetUtil import HandlingError
import Util
import Form
import FormOut
import Progress

def get_hamming_distance(a, b):
    """
    @param a: a sequence of ones and zeros
    @param b: a sequence of ones and zeros
    @return: the hamming distance
    """
    return Util.hamming_distance(a, b)

def get_folded_distance(a, b):
    """
    @param a: a sequence of ones and zeros
    @param b: a sequence of ones and zeros
    @return: a distance
    """
    max_distance = len(a)
    hdist = Util.hamming_distance(a, b)
    return min(hdist, max_distance - hdist)

def get_form():
    """
    @return: the body of a form
    """
    return []

def get_form_out():
    return FormOut.Report()

def get_centers(arrs, ncenters, metric, deadline):
    """
    @param arrs: binary arrays
    @param ncenters: the number of centers to find
    @param metric: this function measures the distance between points
    @param deadline: raise an error if it gets later than this
    @return: the maximum distance to a center, and the list of centers
    """
    assert ncenters > 0
    assert ncenters <= len(arrs)
    # create the distance matrix
    n = len(arrs)
    D = [[metric(a,b) for b in arrs] for a in arrs]
    # get the centers
    best_max_distance = None
    best_centers = None
    for center_indices in itertools.combinations(range(n), ncenters):
        if deadline and time.time() > deadline:
            raise HandlingError('im sorry this calculation has exceeded my attention span')
        max_distance = max(min(D[c][t] for c in center_indices) for t in range(n))
        if not best_centers or max_distance < best_max_distance:
            best_max_distance = max_distance
            best_centers = [arrs[i] for i in center_indices]
    return best_max_distance, best_centers

def get_result_string(arrs, max_distance, centers):
    """
    @param arrs: binary arrays
    @param max_distance: the max dist of a binary string to the nearest center
    @param centers: the centers found by the search
    """
    out = StringIO()
    print >> out, 'number of binary strings:', len(arrs)
    print >> out, 'length of each binary string:', len(arrs[0])
    print >> out, 'maximum distance from a center:', max_distance
    print >> out, 'centers:'
    for center in centers:
        print >> out, ''.join(str(x) for x in center)
    return out.getvalue().strip()

def get_response_content(fs):
    return 'sorry this is not finished'

def get_histogram(arr):
    """
    @param arr: a sequence of numbers
    """
    d = {}
    for value in arr:
        count = d.get(value, 0)
        d[value] = count + 1
    return d

def get_min_squared_error(phenotypes, snp):
    """
    Get the minimum possible squared error for a pair that includes the SNP.
    """
    pass

def main():
    filename = '/home/argriffi/phenotype-challenge/testphen.txt'
    print 'parsing the data'
    # read the input file
    fin = open(filename)
    lines = fin.readlines()
    fin.close()
    # parse the input file
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    rows = [[x.strip() for x in line.split(',')] for line in lines]
    phenotypes = [float(x) for x in rows[0]]
    nphenotypes = len(phenotypes)
    snps = [[int(float(x)) for x in row] for row in rows[1:]]
    for row in snps:
        assert len(row) == nphenotypes
    print 'phenotypes:'
    print phenotypes
    print 'found', nphenotypes, 'phenotypes'
    print 'found', len(snps), 'snps'
    print 'first rows lines of snps:'
    print '\n'.join(''.join(str(x) for x in row) for row in snps[:3])
    # look for uninformative rows and columns
    row_sets = set(tuple(row) for row in snps)
    column_sets = set(tuple(column) for column in zip(*snps))
    print len(row_sets), 'unique snp patterns'
    print len(column_sets), 'unique genetic lines'

def do_break_point_stuff():
    error_list = []
    n = 10
    arr = list(sorted(random.random() for i in range(n)))
    for offset in range(n+1):
        se = 0
        for chunk in (arr[:offset], arr[offset:]):
            if chunk:
                mu = sum(chunk) / len(chunk)
                se += sum((x-mu)**2 for x in chunk)
        error_list.append(se)
    for low, high in zip(error_list[:-1], error_list[1:]):
        print high - low

def do_data_translation():
    """
    Rewrite the data so that it is easily interpretable by C code.
    """
    filename_in = '/home/argriffi/phenotype-challenge/testphen.txt'
    filename_out = '/home/argriffi/phenotype-challenge/testphen-simple.txt'
    # read the input file
    fin = open(filename_in)
    lines = fin.readlines()
    fin.close()
    # parse the input file
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    rows = [[x.strip() for x in line.split(',')] for line in lines]
    phenotypes = [float(x) for x in rows[0]]
    nphenotypes = len(phenotypes)
    snps = [[int(float(x)) for x in row] for row in rows[1:]]
    # validate the data
    for row in snps:
        assert len(row) == nphenotypes
    # center the phenotypes
    mu = sum(phenotypes) / nphenotypes
    phenotypes = [x - mu for x in phenotypes]
    # create the output data
    output_lines = []
    output_lines.extend(str(x) for x in phenotypes)
    output_lines.extend(''.join(str(x) for x in snp) for snp in snps)
    # write the output data
    fout = open(filename_out, 'w')
    fout.write('\n'.join(output_lines) + '\n')
    fout.close()

def do_r_translation():
    """
    Rewrite the data so that it is easily interpretable by R code.
    """
    filename_in = '/home/argriffi/phenotype-challenge/testphen.txt'
    matrix_filename_out = 'matrix.dat'
    vector_filename_out = 'vector.dat'
    # read the input file
    fin = open(filename_in)
    lines = fin.readlines()
    fin.close()
    # parse the input file
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if line]
    rows = [[x.strip() for x in line.split(',')] for line in lines]
    phenotypes = [float(x) for x in rows[0]]
    nphenotypes = len(phenotypes)
    snps = [[int(float(x)) for x in row] for row in rows[1:]]
    # validate the data
    for row in snps:
        assert len(row) == nphenotypes
    # center the phenotypes
    #mu = sum(phenotypes) / nphenotypes
    #phenotypes = [x - mu for x in phenotypes]
    # write the matrix
    flattened_matrix = []
    for row in snps:
        flattened_matrix.extend(row)
    fout = open(matrix_filename_out, 'w')
    print >> fout, ' '.join(str(x) for x in flattened_matrix)
    fout.close()
    # write the vector
    fout = open(vector_filename_out, 'w')
    print >> fout, ' '.join(str(x) for x in phenotypes)
    fout.close()

if __name__ == '__main__':
    do_r_translation()

