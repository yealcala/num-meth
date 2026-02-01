from sage.all import *
from .factorization import *
from .remonte import *
from ...verify import nonNull, matrixMustBeSquare, matrixMustHaveNonNullDiagonal

def inverseDiagonal(D: matrix):
    nonNull(D)
    matrixMustBeSquare(D)
    matrixMustHaveNonNullDiagonal(D)

    n = D.ncols()
    
    D1 = matrix.zero(QQ, n)
    for i in range(n):
        D1[i, i] = 1/D[i, i]
    return D1

def inverseLU(A: matrix):
    """Suponemos A acepta factorización LU"""
    L, U = doolittle(A)
    n = A.ncols()
    A1 = matrix(QQ, n, n)

    for i in range(n):
        ei = vector(QQ, [0]*n)
        ei[i] = 1
        u = remonteArriba(U, remonteAbajo(L, ei))
        A1.set_column(i, u)
    return A1

def solveLU(A: matrix, b: vector) -> vector:
    """Suponemos A acepta factorización LU"""
    L, U = doolittle(A)
    return remonteArriba(U, remonteAbajo(L, b))

def solvePLU(A: matrix, b: vector) -> vector:
    """Suponemos A acepta factorización LU"""
    P, L, U = doolittlePLU(A)
    return remonteArriba(U, remonteAbajo(L, P*b))
