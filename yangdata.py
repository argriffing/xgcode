"""
Read some of the data files from the website of Ziheng Yang.

This does not download the files from the website;
the functions assume that you have the data files available locally.
"""

import numpy
import npcodon

g_mtdna_names = {'human_horai', 'chimp_horai'}

def read_yang_alignments(codons, lines):
    """
    Read the file that Ziheng Yang put on his website.
    For our purposes,
    an alignment will be a list of [taxon_name, codon_index_sequence] pairs.
    @param codons: sequence lowercase codons
    @param lines: text lines of the alignment
    @return: list of alignments
    """
    names = {'hg18', 'pantro2', 'rhemac2', 'mm8', 'rn4'}
    # construct the inverse codon map
    c_to_i = dict((c, i) for i, c in enumerate(codons))
    # read the alignments
    alignments = []
    alignment = []
    expected_len = None
    for line in lines:
        line = line.strip().lower()
        if not line:
            continue
        elements = line.split()
        if elements[0] in names:
            seq = numpy.array([c_to_i[c] for c in elements[1:]], dtype=int)
            if len(seq) * 3 != expected_len:
                raise Exception((len(seq) * 3, expected_len))
            alignment.append((elements[0], seq))
        elif int(elements[0]) == len(names):
            if len(elements) != 2:
                raise Exception(elements)
            expected_len = int(elements[1])
            if alignment:
                alignments.append(alignment)
            alignment = []
        else:
            raise Exception(elements[0])
    if alignment:
        alignments.append(alignment)
    return alignments

def read_yang_mtdna_alignment(codons, lines):
    c_to_i = dict((c, i) for i, c in enumerate(codons))
    name = None
    alignment = []
    segments = []
    for line in lines:
        line = line.strip().lower()
        if line in g_mtdna_names or not line:
            if segments:
                dna = ''.join(segments)
                codons = [dna[i:i+3] for i in range(0, len(dna), 3)]
                seq = numpy.array(
                        [c_to_i[''.join(c)] for c in codons], dtype=int)
                alignment.append((name, seq))
        if line in g_mtdna_names:
            name = line
            segments = []
        elif not line:
            name = None
        else:
            if name:
                segments.append(line)
    return alignment

def get_empirical_summary(ncodons, alignments, t1, t2):
    """
    Get substitution counts.
    It is easy to get empirical codon counts from the substitution counts.
    This is a helper function for get_codon_counts_from_data_files.
    @param ncodons: allowed number of codons
    @param alignments: from the get_yang_alignments function
    @param t1: first taxon name
    @param t2: second taxon name
    @return: subs_counts
    """
    subs_counts = numpy.zeros((ncodons, ncodons), dtype=int)
    for alignment in alignments:
        d = dict(alignment)
        for i, j in zip(d[t1], d[t2]):
            subs_counts[i, j] += 1
    return subs_counts

def get_subs_counts_from_data_files(args):
    """
    @return: numpy ndarray of observed substitution counts
    """
    if args.mtdna:
        code = npcodon.g_code_mito
        stop = npcodon.g_stop_mito
    else:
        code = npcodon.g_code
        stop = npcodon.g_stop
    all_codons = npcodon.enum_codons(stop)
    with open(args.infile) as fin:
        if args.mtdna:
            alignments = [read_yang_mtdna_alignment(all_codons, fin)]
        else:
            alignments = read_yang_alignments(all_codons, fin)
    if args.mtdna:
        t1, t2 = g_mtdna_names
    else:
        t1, t2 = args.t1, args.t2
    subs_counts = get_empirical_summary(64, alignments, t1, t2)
    return subs_counts

