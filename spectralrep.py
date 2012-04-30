"""
This module is about the spectral representation of certain rate matrices.

In particular, it is about the spectral representation of
irreducible reversible rate matrices that represent
continuous time finite state Markov processes.
In this module we never try to find the eigenvalues or eigendecomposition
of asymmetric matrices directly.
Matrices are represented by python numpy arrays.
The scipy linalg eigh and eigvalsh functions are used
because their return values are ordered informatively
and because they use the symmetry constraint to guarantee
eigenvalues with no imaginary component.
"""

import unittest

import numpy as np
import scipy
from scipy import linalg



class TestSpectralRep(unittest.TestCase):

    def test_placeholder(self):
        pass


if __name__ == '__main__':
    unittest.main()

