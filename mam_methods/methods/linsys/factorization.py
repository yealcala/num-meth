from sage.all import *
from ...verify import nonNull, matrixMustBeSquare, matrixMustBeInvertible, matrixMustBeHermitian, matrixMustBeDefinitePositive

def doolittlePLU(A: matrix) -> tuple[matrix, matrix, matrix]:
    nonNull(A)
    matrixMustBeSquare(A)

    n = A.nrows()
    Ac = copy(A)
    P = matrix.identity(QQ, n)
    L = matrix.identity(QQ, n)
    U = matrix.zero(QQ, n)
    
    
    for i in range(0, n):
        piv = i
        for k in range(i + 1, n):
            if abs(Ac[k, i]) > abs(Ac[piv, i]):
                piv = k
        
        if piv != i:
            Ac.swap_rows(i, piv)
            P.swap_rows(i, piv)
            for k in range(i): # Intercambiar en L
                L[i, k], L[piv, k] = L[piv, k], L[i, k]

        for j in range(i, n):
            U[i, j] = A[i, j] - sum([L[i, k]*U[k, j] for k in range(0, i)])
        if U[i,i] != 0:
            for j in range(i+1, n):
                L[j, i] = 1/U[i, i]*(A[j, i] - sum([L[j, k]*U[k, i] for k in range(0, i)]))
    
    return (P, L, U)

def doolittle(A: matrix) -> tuple[matrix, matrix]:
    nonNull(A)
    matrixMustBeSquare(A)
    matrixMustBeInvertible(A)

    L = matrix.identity(QQ, A.nrows())
    U = matrix.zero(QQ, A.nrows())
    n = A.nrows()
    
    for i in range(0, n):
        for j in range(i, n):
            U[i, j] = A[i, j] - sum([L[i, k]*U[k, j] for k in range(0, i)])
        if U[i,i] == 0:
            raise ZeroDivisionError(f"No se puede factorizar: menor principal nulo en índice {i}.")
        for j in range(i+1, n):
            L[j, i] = 1/U[i, i]*(A[j, i] - sum([L[j, k]*U[k, i] for k in range(0, i)]))
    
    return (L, U)

def crout(A: matrix) -> tuple[matrix, matrix]:
    nonNull(A)
    matrixMustBeSquare(A)
    matrixMustBeInvertible(A)

    L = matrix.zero(QQ, A.nrows())
    U = matrix.identity(QQ, A.nrows())
    n = A.nrows()

    for i in range(0, n):
        for j in range(i, n):
            L[j, i] = A[j, i] - sum(L[j, k] * U[k, i] for k in range(i))
        
        if L[i, i] == 0:
            raise ZeroDivisionError(f"No se puede factorizar: menor principal nulo en índice {i}.")

        for j in range(i + 1, n):
            U[i, j] = (1 / L[i, i]) * (A[i, j] - sum(L[i, k] * U[k, j] for k in range(i)))
    return (L, U)

def cholesky(A: matrix) -> tuple[matrix, matrix]:
    nonNull(A)
    matrixMustBeSquare(A)
    matrixMustBeInvertible(A)
    matrixMustBeHermitian(A)
    matrixMustBeDefinitePositive(A)

    n = A.nrows()
    B = matrix(SR, n, n)
    
    for k in range(n):
        B[k, k] = sqrt(A[k, k] - sum([B[k, r]**2 for r in range(k)]))
        for i in range(k+1, A.nrows()):
            B[i, k] = 1/B[k, k]*(A[i, k] - sum([B[i, r]*B[k, r] for r in range(k)]))

    
    return (B, B.transpose())