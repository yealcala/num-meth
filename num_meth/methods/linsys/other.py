# num_meth/methods/linsys/other.py

"""Provides other methods.

The module contains the following functions:

- `inverseDiagonal(D)` - Calculates the inverse of a diagonal matrix.
- `inverseLU(A)` - Calculates the inverse of a matrix using LU factorization.
- `solveLU(A, b)` - Solves Au = b linear system using LU factorization.
- `solvePLU(A, b)` - Solves Au = b linear system using PA = LU factorization.

"""
from sage.all import matrix, vector
from .factorization import *
from .remonte import *
from ...verify import nonNull, matrixMustBeSquare, matrixMustHaveNonNullDiagonal

def inverseDiagonal(D: matrix):
    """Calculates the inverse of a diagonal matrix.        

    Args:
        D (matrix): A square diagonal n-matrix with all diagonal elements non null.

    Returns:
        matrix: The inverse of D.

    Raises:
        ValueError: If some argument is not valid.
    """
    nonNull(D)
    matrixMustBeSquare(D)
    matrixMustHaveNonNullDiagonal(D)

    n = D.ncols()
    
    D1 = matrix.zero(QQ, n)
    for i in range(n):
        D1[i, i] = 1/D[i, i]
    return D1

def inverseLU(A: matrix):
    """Calculates the inverse of a matrix using A = LU factorization.        

    Args:
        A (matrix): A square diagonal n-matrix which accepts unique LU factorization.

    Returns:
        matrix: The inverse of A.

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If the linear system has multiple or none solutions.
    """    
    L, U = doolittle(A)
    n = A.ncols()
    A1 = matrix(QQ, n, n)

    for i in range(n):
        ei = vector(QQ, [0]*n)
        ei[i] = 1
        u = upperBackSubs(U, lowerBackSubs(L, ei))
        A1.set_column(i, u)
    return A1

def solveLU(A: matrix, b: vector) -> vector:
    """Calculates the solution of Au = b linear matrix using A = LU factorization.        

    Args:
        A (matrix): A square diagonal n-matrix which accepts unique LU factorization.
        b (vector)): A n-vector.

    Returns:
        vector: The solution of Au = b linear system.

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If the linear system has multiple or none solutions.
    """  
    L, U = doolittle(A)
    return upperBackSubs(U, lowerBackSubs(L, b))

def solvePLU(A: matrix, b: vector) -> vector:
    """Calculates the solution of Au = b linear matrix using PA = LU factorization.        

    Args:
        A (matrix): A square diagonal n-matrix which accepts unique P-LU factorization.
        b (vector)): A n-vector.

    Returns:
        vector: The solution of Au = b linear system.

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If the linear system has multiple or none solutions.
    """  
    P, L, U = doolittlePLU(A)
    return upperBackSubs(U, lowerBackSubs(L, P*b))
