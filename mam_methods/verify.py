from sage.all import *
from sage.structure.element import Matrix, Vector

def nonNull(*args):
    for i, arg in enumerate(args):
        if arg is None:
            raise ValueError(f"None argument found in position: {i}")
    return True

# <----------- MATRIX ----------->

def objMustBeMatrix(o: any) -> None:
    if not isinstance(o, Matrix):
        raise ValueError("Object must be Matrix!")

def matrixMustBeSquare(A: matrix) -> None:
    objMustBeMatrix(A)
    nonNull(A)
    if A.ncols() != A.nrows():
        raise ValueError("Matrix must be squared!")
    
def matrixMustBeInvertible(A: matrix) -> None:
    matrixMustBeSquare(A)
    if not A.is_invertible():
        raise ValueError("Matrix must be invertible!")

def matrixMustBeSymmetric(A: matrix) -> None:
    matrixMustBeSquare(A)
    if A != A.transpose():
        raise ValueError("Matrix must be symmetric!")

def matrixMustBeHermitian(A: matrix) -> None:
    matrixMustBeSquare(A)
    if A != A.conjugate_transpose():
        raise ValueError("Matrix must be hermitian!")
    
def matrixMustBeDefinitePositive(A: matrix) -> None:
    matrixMustBeSquare(A)
    if A != A.conjugate_transpose():
        raise ValueError("Matrix must be definite positive!")
    
def matrixMustHaveNonNullDiagonal(A: matrix) -> None:
    matrixMustBeSquare(A)
    for i in range(0, A.ncols()):
        if(A[i,i] == 0):
            raise ValueError("Matrix must have non-null diagonal!")
        
def matrixMustBeTriangular(A: matrix, upper: bool = False) -> None:
    matrixMustBeSquare(A)
    side = "upper" if upper else "lower"
    if not A.is_triangular(side):
        raise ValueError(f"A must be {side} triangular!")

# <----------- VECTOR ----------->
def objMustBeVector(o: any) -> None:
    if not isinstance(o, Vector):
        raise ValueError("Object must be Vector!")
    
def vectorMustBeDimension(v: vector, n: int) -> None:
    objMustBeVector(v)
    nonNull(v, n)
    if len(list(v)) != n:
        raise ValueError(f"Vector must be of dimension {n}")
    