"""
Try to do some naive symbolic matrix stuff.
"""

from StringIO import StringIO
import numpy as np
import scipy
import scipy.linalg

import Form
import FormOut


#TODO actually use precedence to place parentheses
# precedence:
# 0: negation, inversion, pseudoinversion
# 1: product
# 2: sum

# Terminals:
# zero
# one
# column vector of ones with symbolic integer length
# row vector of ones with symbolic integer length
# identity matrix with symbolic integer size
# symbolic integer
# symbolic real scalar
# symbolic column vector of reals with symbolic integer length
# symbolic row vector of reals with symbolic integer length
# symbolic symmetric matrix of reals with symbolic integer size

# everything has a two dimensional shape

class Terminal: pass

class Zero(Terminal):
    def latex(self):
        return '0'
    def get_shape(self):
        return (One(), One())

class One(Terminal):
    def latex(self):
        return '1'
    def get_shape(self):
        return (One(), One())

class OnesColumnVector(Terminal):
    def __init__(self, length):
        """
        @param length: SymbolicInteger length
        """
        self.length = length
    def latex(self):
        return 'e_{%s}' % self.length.latex()
    def get_shape(self):
        return (self.length, One())

class OnesRowVector(Terminal):
    def __init__(self, length):
        """
        @param length: SymbolicInteger length
        """
        self.length = length
    def latex(self):
        return 'e^T_{%s}' % self.length.latex()
    def get_shape(self):
        return (One(), self.length)

class Identity(Terminal):
    def __init__(self, order):
        """
        @param order: SymbolicInteger order
        """
        self.order = order
    def latex(self):
        return 'I_{%s}' % self.order.latex()
    def get_shape(self):
        return (self.order, self.order)

class SymbolicInteger(Terminal):
    def __init__(self, name):
        self.name = name
    def latex(self):
        return self.name
    def get_shape(self):
        return (One(), One())

class SymbolicScalar(Terminal):
    def __init__(self, name):
        self.name = name
    def latex(self):
        return self.name
    def get_shape(self):
        return (One(), One())

class SymbolicColumnVector(Terminal):
    def __init__(self, name, length):
        """
        @param name: string name
        @param length: SymbolicInteger length
        """
        self.name = name
        self.length = length
    def latex(self):
        return self.name
    def get_shape(self):
        return (self.length, One())

class SymbolicRowVector(Terminal):
    def __init__(self, name, length):
        """
        @param name: string name
        @param length: SymbolicInteger length
        """
        self.name = name
        self.length = length
    def latex(self):
        return '%s^T' % self.name
    def get_shape(self):
        return (One(), self.length)

class SymbolicSymmetricMatrix(Terminal):
    def __init__(self, name, order):
        """
        @param name: string name
        @param order: SymbolicInteger order
        """
        self.name = name
        self.order = order
    def latex(self):
        return self.name
    def get_shape(self):
        return (self.order, self.order)


# define non-terminal things

class UnaryOperator:
    def get_principal_shape(self):
        return self.element.get_principal_shape()
    def get_block_structure(self):
        return self.element.get_block_structure()

class Inverse(UnaryOperator):
    """
    Multiplicative inverse of a single thing.
    Inverse of integer or scalar or matrix or sum or product or inverse.
    """
    precedence = 0
    def __init__(self, element):
        self.element = element
    def latex(self):
        #return '\\dfrac{1}{%s}' % self.element.latex()
        return '{%s}^{-1}' % self.element.latex()
    def get_shape(self):
        return self.element.get_shape()

class Negative(UnaryOperator):
    """
    Negative of a single thing.
    """
    precedence = 0
    def __init__(self, element):
        self.element = element
    def latex(self):
        return '-%s' % self.element.latex()
    def expand(self):
        try:
            self.element.distribute_negation()
        except AttributeError, e:
            pass
    def get_shape(self):
        return self.element.get_shape()

class InverseInH(UnaryOperator):
    """
    This is a specific pseudoinverse.
    It assumes that an order n symmetric matrix is rank n-1
    and that its eigenvector with zero eigenvalue is the constant vector.
    That is, it assumes the matrix is doubly centered.
    """
    precedence = 0
    def __init__(self, element):
        self.element = element
    def latex(self):
        return '{%s}^\ddag' % self.element.latex()
    def get_shape(self):
        return self.element.get_shape()

class PinvProj(UnaryOperator):
    """
    Pseudoinverse of double centering.
    The element is a matrix or block matrix or higher order thing.
    """
    precedence = 0
    def __init__(self, element):
        self.element = element
    def latex(self):
        return '\\left( H %s H \\right)^\ddag' % self.element.latex()
    def expanded(self):
        block_structure = None
        try:
            block_structure = self.element.get_block_structure()
        except AttributeError, e:
            pass
        if block_structure:
            # use the block structure to make block matrices
            H = get_block_H_expression(block_structure)
            P = get_block_P_expression(block_structure)
        else:
            # no block structure is available
            order = self.get_shape()[0]
            H = get_H_expression(order)
            P = get_P_expression(order)
        terma = Inverse(Sum([Product([H, self.element, H]), P]))
        termb = Negative(P)
        return Sum([terma, termb])
    def get_shape(self):
        return self.element.get_shape()

class Product:
    """
    Product of multiple things.
    """
    precedence = 1
    def __init__(self, elements):
        self.elements = elements
    def latex(self):
        inside = ' '.join(x.latex() for x in self.elements)
        return '\\left( %s \\right)' % inside
    def get_shape(self):
        #TODO this is subtle and probably wrong
        first_shape = self.elements[0]
        last_shape = self.elements[-1]
        return (first_shape[0], last_shape[-1])
    def get_principal_shape(self):
        #TODO this is subtle and probably wrong
        try:
            return self.elements[0].get_principal_shape()
        except AttributeError, e:
            return self.get_shape()
    def get_block_structure(self):
        #TODO this is subtle and probably wrong
        try:
            return self.elements[0].get_block_structure()
        except AttributeError, e:
            return self.get_shape()

class Schur:
    """
    Schur complement of lower right block in a 2x2 block matrix.
    The single member element is the 2x2 block matrix.
    """
    precedence = 1
    def __init__(self, element):
        self.element = element
    def latex(self):
        return '%s / \\cdot' % self.element.latex()
    def get_shape(self):
        return self.element.get_principal_shape()
    def expanded(self):
        """
        Assume the element is a two by two block matrix.
        """
        #TODO do not make assumptions about the member
        (a, b), (c, d) = self.element.elements
        return Sum([a, Negative(Product([b, Inverse(c), d]))])

class Sum:
    """
    Sum of multiple things.
    """
    precedence = 2
    def __init__(self, elements):
        self.elements = elements
    def latex(self):
        parts = []
        for i, x in enumerate(self.elements):
            if isinstance(x, Negative):
                parts.append('-')
                parts.append(x.element.latex())
            else:
                if i:
                    parts.append('+')
                parts.append(x.latex())
        return '\\left( %s \\right)' % ' '.join(parts)
    def distribute_negation(self):
        next_elements = []
        for e in self.elements:
            if isinstance(e, Negative):
                next_elements.append(e.element)
            else:
                next_elements.append(Negative(e))
    def get_shape(self):
        return self.elements[0].get_shape()
    def get_principal_shape(self):
        return self.elements[0].get_principal_shape()
    def get_block_structure(self):
        return self.elements[0].get_block_structure()


# define non-terminal block matrices of non-symbolic fixed size

class BlockMatrixTwoByTwo:
    """
    Two by two block matrix.
    Yes this is quite hardcoded.
    Assume that the block structure is symmetric,
    so that it can be defined by two numbers.
    Elements are in a row major list of lists.
    """
    def __init__(self, elements):
        self.elements = elements
    def latex(self):
        (a, b), (c, d) = self.elements
        out = StringIO()
        print >> out, '\\begin{pmatrix}'
        print >> out, a.latex(), '&', b.latex(), '\\\\'
        print >> out, c.latex(), '&', d.latex()
        print >> out, '\\end{pmatrix}'
        return out.getvalue()
    def get_shape(self):
        (a, b), (c, d) = self.elements
        nrows = Sum([a.get_shape()[0], c.get_shape()[0]])
        ncols = Sum([a.get_shape()[1], b.get_shape()[1]])
        return (nrows, ncols)
    def get_principal_shape(self):
        (a, b), (c, d) = self.elements
        return a.get_shape()
    def get_block_structure(self):
        """
        Assume that the block structure can be defined by two numbers.
        """
        (A, B), (C, D) = self.elements
        # Assume that block A is square and so therefore is block D.
        return (A.get_shape()[0], D.get_shape()[0])


# helper functions

def get_H_expression(order):
    """
    H is a centering matrix.
    """
    P = get_P_expression(order)
    return Sum([Identity(order), Negative(P)])

def get_P_expression(order):
    """
    P is (I - H) where H is the centering matrix.
    """
    return Product([
        OnesColumnVector(order),
        Inverse(Product([OnesRowVector(order), OnesColumnVector(order)])),
        OnesRowVector(order)])

def get_block_H_expression(block_structure):
    """
    H is a centering matrix.
    """
    #TODO automatically expand this from the simpler H expression
    na, nb = block_structure
    prefix = Negative(Inverse(Sum([na, One()])))
    A = Sum([
        Identity(na), 
        Product([prefix, OnesColumnVector(na), OnesRowVector(na)])])
    B = Product([prefix, OnesColumnVector(na)])
    C = Product([prefix, OnesRowVector(na)])
    D = Sum([One(), prefix])
    elements = [[A, B], [C, D]]
    return BlockMatrixTwoByTwo(elements)

def get_block_P_expression(block_structure):
    """
    P is (I - H) where H is the centering matrix.
    """
    #TODO automatically expand this from the simpler P expression
    na, nb = block_structure
    prefix = Inverse(Sum([na, One()]))
    A = Product([prefix, OnesColumnVector(na), OnesRowVector(na)])
    B = Product([prefix, OnesColumnVector(na)])
    C = Product([prefix, OnesRowVector(na)])
    D = prefix
    elements = [[A, B], [C, D]]
    return BlockMatrixTwoByTwo(elements)


def get_form():
    """
    @return: the body of a form
    """
    # define the form objects
    form_objects = []
    return form_objects

def get_form_out():
    return FormOut.Report()

def get_response_content(fs):
    # define lhs and rhs math objects
    n = SymbolicInteger('n')
    lhs = PinvProj(SymbolicSymmetricMatrix('A', n)).expanded()
    rhs = Schur(PinvProj(
        BlockMatrixTwoByTwo([
            [
                SymbolicSymmetricMatrix('A', n),
                SymbolicColumnVector('b', n)],
            [
                SymbolicRowVector('b', n),
                SymbolicScalar('c')]])).expanded())
    # get some tex code
    out = StringIO()
    print >> out, '\\begin{equation*}'
    print >> out, '\\text{lhs} ='
    print >> out, lhs.latex()
    print >> out, '\\end{equation*}'
    print >> out
    print >> out, '\\begin{equation*}'
    print >> out, '\\text{rhs} ='
    print >> out, rhs.latex()
    print >> out, '\\end{equation*}'
    return out.getvalue()

