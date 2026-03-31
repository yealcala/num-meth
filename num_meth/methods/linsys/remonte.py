# num_meth/methods/linsys/remonte.py

"""Provides the back substitution method.

The module contains the following functions:

- `upperBackSubs(A, b, tol)` - Applies upper back substitution to linear system Au = b with A upper triangular.
- `lowerBackSubs(A, b, tol)` - Applies lower back substitution to linear system Au = b with A lower triangular.

"""

from sage.all import matrix, vector, QQ, oo
from ...exceptions import ResidualToleranceError
from ...utils import pnorm
from ...verify import nonNull, matrixMustHaveNonNullDiagonal, matrixMustBeTriangular, vectorMustBeDimension

def __checkSol(A: matrix, b: vector, u: vector, tol: float) -> bool:
    res = A * u - b
    if all(abs(r) < tol for r in res):
        print("La solución es consistente.")
        return True
    else:
        print("El residuo es significativo.")
        raise ResidualToleranceError(pnorm(res, oo))
    
def upperBackSubs(A: matrix, b: vector, tol: float | None = None) -> vector:

    """Resolves Au = b linear system applying upper back substitution.

    Examples:
        >>> A = matrix(QQ, 3, [1, 2, 3, 0, 2, 3, 0, 0, 3])
        >>> b = vector(QQ, [1, 2, 3])

        >>> show(remonteArriba(A, b))
        

    Args:
        A (matrix): A square n-matrix being upper triangular.
        b (vector): A n-vector
        tol (Optional[float]): Maximum error that given solution can have

    Returns:
        vector: The solution of the linear system Au = b.

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If the linear system has multiple or none solutions.
        ResidualToleranceError: If given solution's and real solution's difference is greater than tolerance.
    """
    nonNull(A, b)
    matrixMustBeTriangular(A, upper=True)
    matrixMustHaveNonNullDiagonal(A)

    n = A.nrows()
    vectorMustBeDimension(b, n)

    u = vector(QQ, [0]*n) 

    for i in range(n-1, -1, -1):
        if A[i, i] == 0:
            raise ZeroDivisionError(f"El sistema no tiene solución única (A[{i},{i}] es nulo).")
        u[i] = 1/A[i,i]*(b[i] - sum([(u[k]*A[i,k]) for k in range(i+1, n)]))
        
    if tol is not None: __checkSol(A, b, u, tol)
    return u



def lowerBackSubs(A: matrix, b: vector, tol: float | None = None) -> vector:

    """Resolves Au = b linear system applying lower back substitution.

    Examples:
        >>> A = matrix(QQ, 3, [1, 2, 3, 0, 2, 3, 0, 0, 3]).transpose()
        >>> b = vector(QQ, [1, 2, 3])

        >>> show(remonteArriba(A, b))
        

    Args:
        A (matrix): A square n-matrix being lower triangular.
        b (vector): A n-vector
        tol (Optional[float]): Maximum error that given solution can have

    Returns:
        vector: The solution of the linear system Au = b.

    Raises:
        ValueError: If some argument is not valid.
        ZeroDivisionError: If the linear system has multiple or none solutions.
        ResidualToleranceError: If given solution's and real solution's difference is greater than tolerance.
    """
    nonNull(A, b)
    matrixMustBeTriangular(A, upper=False)
    matrixMustHaveNonNullDiagonal(A)
    
    n = A.nrows()
    vectorMustBeDimension(b, n)

    u = vector(QQ, [0]*n) 
    
    for i in range(n):
        if A[i, i] == 0:
            raise ZeroDivisionError(f"El sistema no tiene solución única (A[{i},{i}] es nulo).")
        u[i] = 1/A[i,i]*(b[i] - sum([(u[k]*A[i,k]) for k in range(i)]))
    
    if tol is not None: __checkSol(A, b, u, tol)
    return u