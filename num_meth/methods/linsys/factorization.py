# num_meth/methods/linsys/calculations.py

"""Provides methods for factorizing a matrix.

The module contains the following functions:

- `doolittlePLU(A)` - Returns matrices for PA = LU factorization where all elements in L's diagonal are 1.
- `doolittle(A)` - Returns matrices for A = LU factorization where all elements in L's diagonal are 1.
- `crout(A)` - Returns matrices for A = LU factorization where all elements in U's diagonal are 1.
- `cholesky(A)` - Returns matrices for A = BB^t factorization where B is lower triangular.
"""

from sage.all import matrix, copy, QQ, SR, sqrt
from ...verify import nonNull, matrixMustBeSquare, matrixMustBeInvertible, matrixMustBeHermitian, matrixMustBeDefinitePositive

def doolittlePLU(A: matrix) -> tuple[matrix, matrix, matrix]:
    """Calculates PA = LU factorization of Doolittle (L's diagonal's elements are 1)

    Examples:
        >>> A = matrix(QQ, [[1, -1, 1, -1], [-1, 3, -3, 3], [2, -4, 7, -7], [-3, 7, -10, 14]])
        >>> show("Doolittle, PA = LU:", doolittlePLU(A))
        

    Args:
        A (matrix): A square n-matrix with all leading principal minors non-null

    Returns:
        tuple: A tuple with three components,
            - 0 (matrix): The permutation matrix P
            - 1 (matrix): The matrix L with all it's diagonal elements being 1
            - 2 (matrix): The matrix U

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If there is a non-null leading principal minor.
    """
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

    """Calculates A = LU factorization of Doolittle (L's diagonal's elements are 1)

    Examples:
        >>> A = matrix(QQ, [[1, -1, 1, -1], [-1, 3, -3, 3], [2, -4, 7, -7], [-3, 7, -10, 14]])
        >>> show("Doolittle, A = LU:", doolittle(A))
        

    Args:
        A (matrix): A square n-matrix with all leading principal minors non-null

    Returns:
        tuple: A tuple with two components,
            - 0 (matrix): The matrix L with all it's diagonal elements being 1
            - 1 (matrix): The matrix U

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If there is a non-null leading principal minor.
    """
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
    """Calculates A = LU factorization of Doolittle (U's diagonal's elements are 1)

    Examples:
        >>> A = matrix(QQ, [[1, -1, 1, -1], [-1, 3, -3, 3], [2, -4, 7, -7], [-3, 7, -10, 14]])
        >>> show("Crout, A = LU:", crout(A))
        

    Args:
        A (matrix): A square n-matrix with all leading principal minors non-null

    Returns:
        tuple: A tuple with two components,
            - 0 (matrix): The matrix L
            - 1 (matrix): The matrix U with all it's diagonal elements being 1

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If there is a non-null leading principal minor.
    """
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
    """Calculates A = BB^t factorization of Cholesky

    Examples:
        >>> A = matrix(QQ, 3, [4, 2, 1, 2, 5, -2, 1, -2, 7])
        >>> show("Cholesky:", cholesky(A))
        

    Args:
        A (matrix): A hermitian invertible square n-matrix being definite positive.

    Returns:
        tuple: A tuple with two components,
            - 0 (matrix): The matrix B being lower triangular.
            - 1 (matrix): The matrix B^t being upper triangular and the transpose of B.

    Raises:
        ValueError: If some argument is not valid.
    """
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