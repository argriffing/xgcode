"""Compute the algebraic connectivity of a Schur complement tree.
"""

from StringIO import StringIO

import scipy

import Form
import FormOut
import Ftree
import FtreeIO

g_default_tree = '(((3:0.333333333333, 4:0.5)7:1.0, 5:1.0)8:0.088383145868, (1:1.0, 2:0.5)6:0.911616854132)r;'

def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = [
            Form.MultiLine('tree', 'tree', g_default_tree)]
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    T, B, N = FtreeIO.newick_to_TBN(fs.tree)
    leaves = Ftree.T_to_leaves(T)
    L_schur = Ftree.TB_to_L_schur(T, B, leaves)
    mu = scipy.linalg.eigh(L_schur, eigvals_only=True)[1]
    return str(mu)
