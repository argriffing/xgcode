"""For now define the fitch66 distances between amino acids.
"""

import unittest
import StringIO

g_fitch66_distance = {}

# I think this gives the parsimonious number of changes between the amino acids.
g_fitch66_multiline_data_string = """
H FITW660101
D Mutation values for the interconversion of amino acid pairs (Fitch, 1966)
R PMID:5917736
A Fitch, W.M.
T An improved method of testing for evolutionary homology
J J. Mol. Biol. 16, 9-16 (1966)
M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV
      0.
      2.      0.
      2.      2.      0.
      1.      2.      1.      0.
      2.      1.      2.      2.      0.
      2.      1.      2.      2.      2.      0.
      1.      2.      2.      1.      3.      1.      0.
      1.      1.      2.      1.      1.      2.      1.      0.
      2.      1.      1.      1.      2.      1.      2.      2.      0.
      2.      2.      1.      2.      2.      3.      3.      2.      2.      0.
      2.      1.      2.      2.      2.      1.      2.      1.      1.      1.      0.
      2.      1.      1.      2.      3.      1.      1.      2.      2.      2.      2.      0.
      2.      1.      2.      3.      3.      2.      2.      2.      3.      1.      1.      1.      0.
      2.      2.      2.      2.      1.      2.      3.      2.      2.      1.      1.      3.      2.      0.
      1.      1.      2.      2.      2.      1.      2.      2.      1.      2.      1.      2.      2.      2.      0.
      1.      1.      1.      2.      1.      1.      2.      1.      2.      1.      1.      2.      2.      1.      1.      0.
      1.      1.      1.      2.      2.      2.      2.      2.      2.      1.      2.      1.      1.      2.      1.      1.      0.
      2.      1.      3.      3.      1.      1.      2.      1.      3.      3.      1.      2.      2.      2.      2.      1.      2.      0.
      2.      2.      1.      1.      1.      1.      2.      2.      1.      2.      2.      2.      3.      1.      2.      1.      2.      2.      0.
      1.      2.      2.      1.      1.      2.      1.      1.      2.      1.      1.      2.      1.      1.      2.      2.      2.      2.      2.      0.
"""

def _init_fitch66_distance():
    """
    Initialize the global dictionary that gives distances between amino acids.
    This function basically just reads the fitch66 multiline data string.
    """
    # this is the hardcoded ordering of the amino acids in the matrix lines
    ordered_letters = list('ARNDCQEGHILKMFPSTWYV')
    # declare that we want to modify this global variable
    global g_fitch66_distance
    # get processed lines
    raw_lines = StringIO.StringIO(g_fitch66_multiline_data_string).readlines()
    stripped_lines = [line.strip() for line in raw_lines]
    nonempty_stripped_lines = [line for line in stripped_lines if line]
    # read the matrix from the processed lines
    header_line_count = 7
    matrix_lines = nonempty_stripped_lines[header_line_count:]
    assert len(matrix_lines) == len(ordered_letters)
    for first_aa_letter, matrix_line in zip(ordered_letters, matrix_lines):
        distances = [float(v) for v in matrix_line.split()]
        for second_aa_index, distance in enumerate(distances):
            second_aa_letter = ordered_letters[second_aa_index]
            g_fitch66_distance[(first_aa_letter, second_aa_letter)] = distance
            g_fitch66_distance[(second_aa_letter, first_aa_letter)] = distance

def get_fitch66_distance(first_aa, second_aa):
    """
    This function wraps the global dictionary.
    @param first_aa: a capitalized amino acid letter
    @param second_aa: another capitalized amino acid letter
    @return: the fitch66 distance between the amino acids
    """
    return g_fitch66_distance[(first_aa, second_aa)]

# initialize the fitch66 distance dictionary
_init_fitch66_distance()


class TestAminoAcid(unittest.TestCase):

    def test_fitch66(self):
        """
        Do a sanity check of the fitch66 distance table.
        """
        # assert that the distance dictionary is symmetric
        for first_aa in 'ARNDCQEGHILKMFPSTWYV':
            for second_aa in 'ARNDCQEGHILKMFPSTWYV':
                first_distance = g_fitch66_distance[(first_aa, second_aa)]
                second_distance = g_fitch66_distance[(second_aa, first_aa)]
                self.assertEqual(first_distance, second_distance)
        # assert that there is no other cheese in the distance dictionary
        self.assertEqual(len(g_fitch66_distance), 20*20)
        # check a couple of entries
        self.assertEqual(get_fitch66_distance('A', 'N'), 2.0)
        self.assertEqual(get_fitch66_distance('D', 'A'), 1.0)


def main():
    """
    Run tests by default.
    """
    suite = unittest.TestLoader().loadTestsFromTestCase(TestAminoAcid)
    unittest.TextTestRunner(verbosity=2).run(suite)

if __name__ == '__main__':
    main()

