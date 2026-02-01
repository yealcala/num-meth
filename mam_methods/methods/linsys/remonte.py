from sage.all import *
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
    
def remonteArriba(A: matrix, b: vector, tol: float | None = None) -> vector:
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



def remonteAbajo(A: matrix, b: vector, tol: float | None = None) -> vector:
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