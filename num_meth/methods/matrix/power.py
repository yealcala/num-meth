# num_meth/methods/matrix/power.py

"""Provides methods for power iteration of matrices.

The module contains the following functions:

- `powerIteration(A, z0, tol, itMax, norm)` - Power iteration.
- `shiftedPowerIteration(A, d, z0, tol, itMax, norm)` - Shifted power iteration.
- `inversePowerIteration(A, z0, tol, itMax, norm)` - Power iteration.
- `shiftedInversePowerIteration(A, d, z0, tol, itMax, norm)` - Shifted power iteration.
- `rayleigh(A, x0, tol, itMax)` - Shifted power iteration.
"""
from sage.all import *
from ...utils import pnorm
from ..linsys.factorization import doolittlePLU
from ..linsys.remonte import upperBackSubs, lowerBackSubs
from ...verify import nonNull, matrixMustBeSquare, matrixMustBeSymmetric, vectorMustBeDimension


def powerIteration(A: matrix, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo) -> tuple[float, vector, int]:
    """Power iteration.

    Examples:
        >>> A = matrix(QQ, 3, [3, 2, 1, 2, 3, 0, 1, 0, 3])
        >>> x0 = vector(QQ, [1, 1, 1])
        >>> show("Max:", powerIteration(A, x0, tol=1e-10)[0].n())
        

    Args:
        A (matrix): An square n-matrix.
        z0 (vector): Initial vector for iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        itMax (Optional[int]): Maximum iterations allowed. Default: 100.
        norm (Optional[int]): Norm to be used. Default: oo.

    Returns:
        tuple:
            - 0: Approximated maximum absolute eigenvalue.
            - 1: Approximated eigenvector related to the eigenvalue.
            - 2: Number of iterations made.
    """
    nonNull(A, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())

    z = z0.change_ring(A.base_ring())/pnorm(z0, norm)

    lOld = 0
    k = 0
    error = tol + 1
    while error > tol and k < itMax:
        y = A*z
        zOld = z

        normaY = pnorm(y, norm)
        if normaY == 0: return (0, z, k)
        z = y / normaY

        j = 0
        maxVal = -1
        for i in range(len(z)):
            if abs(z[i]) > maxVal:
                maxVal = abs(z[i])
                j = i

        lCurr = y[j] / zOld[j]

        error = abs(lCurr - lOld)
        lOld = lCurr
        k += 1

    return (lOld, z, k)

def shiftedPowerIteration(A: matrix, d: float, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo):
    """Shifted power iteration.

    Examples:
        >>> A = matrix(QQ, 3, [3, 2, 1, 2, 3, 0, 1, 0, 3])
        >>> x0 = vector(QQ, [1, 1, 1])
        >>> show("Furthest to 3:", shiftedPowerIteration(A, 3, x0, tol=1e-10)[0].n())
        

    Args:
        A (matrix): An square n-matrix.
        d (float): A number.
        z0 (vector): Initial vector for iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        itMax (Optional[int]): Maximum iterations allowed. Default: 100.
        norm (Optional[int]): Norm to be used. Default: oo.

    Returns:
        tuple:
            - 0: Approximated eigenvalue furthest from number d in absolute.
            - 1: Approximated eigenvector related to the eigenvalue.
            - 2: Number of iterations made.
    """
    nonNull(A, d, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())
    n = A.nrows()

    As = A - d * matrix.identity(A.base_ring(), n)
    mu, vep, k = powerIteration(As, z0, tol, itMax, norm)
    
    return (mu + d, vep, k)

def inversePowerIteration(A: matrix, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo):
    """Inverse power iteration.

    Examples:
        >>> A = matrix(QQ, 3, [3, 2, 1, 2, 3, 0, 1, 0, 3])
        >>> x0 = vector(QQ, [1, 1, 1])
        >>> show("Min:", inversePowerIteration(A, x0, tol=1e-10)[0].n())


    Args:
        A (matrix): An square n-matrix.
        z0 (vector): Initial vector for iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        itMax (Optional[int]): Maximum iterations allowed. Default: 100.
        norm (Optional[int]): Norm to be used. Default: oo.

    Returns:
        tuple:
            - 0: Approximated minimum absolute eigenvalue.
            - 1: Approximated eigenvector related to the eigenvalue.
            - 2: Number of iterations made.
    """
    nonNull(A, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())
    z = z0.change_ring(A.base_ring()) / pnorm(z0, norm)

    P, L, U = doolittlePLU(A) 
    
    lOld = 0
    k = 0
    error = tol + 1
    while error > tol and k < itMax:
        y = upperBackSubs(U, lowerBackSubs(L, P*z))
        zOld = z
        
        normaY = pnorm(y, norm)
        if normaY == 0: return (0, z, k)
        z = y / normaY

        j = 0
        maxVal = -1
        for i in range(len(z)):
            if abs(z[i]) > maxVal:
                maxVal = abs(z[i])
                j = i

        lCurr =  zOld[j] / y[j]

        error = abs(lCurr - lOld)
        lOld = lCurr
        k += 1

    return (lOld, z, k)

def shiftedInversePowerIteration(A: matrix, d: float, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo):
    """Shifted inverse power iteration.

    Examples:
        >>> A = matrix(QQ, 3, [3, 2, 1, 2, 3, 0, 1, 0, 3])
        >>> x0 = vector(QQ, [1, 1, 1])
        >>> show("Nearest to 2:", shiftedInversePowerIteration(A, 2, x0, tol=1e-10, itMax=500)[0].n())
        

    Args:
        A (matrix): An square n-matrix.
        d (float): A number.
        z0 (vector): Initial vector for iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        itMax (Optional[int]): Maximum iterations allowed. Default: 100.
        norm (Optional[int]): Norm to be used. Default: oo.

    Returns:
        tuple:
            - 0: Approximated eigenvalue nearest to number d in absolute.
            - 1: Approximated eigenvector related to the eigenvalue.
            - 2: Number of iterations made.
    """
    nonNull(A, d, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())
    n = A.nrows()

    As = A - d * matrix.identity(A.base_ring(), n)
    mu, vep, k = inversePowerIteration(As, z0, tol, itMax, norm)
    
    return (mu + d, vep, k)

def rayleigh(A: matrix, x0, tol: float = 1e-10, itMax: int = 100):
    """Rayleigh

    Examples:
        >>> A = matrix(RR, [[10,  2,  1], [ 2, 15, -2], [ 1, -2, 20]])
        >>> x0 = vector(RR, [1, 1, 1])
        >>> show(rayleigh(A, x0, tol=1e-12))
        

    Args:
        A (matrix): A symmetrical square n-matrix.
        x0 (vector): Initial vector for iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        itMax (Optional[int]): Maximum iterations allowed. Default: 100.

    Returns:
        tuple:
            - 0: Approximated eigenvalue.
            - 1: Approximated eigenvector related to the eigenvalue.
            - 2: Number of iterations made.
    """
    nonNull(A, x0, tol, itMax)
    matrixMustBeSymmetric(A)
    
    x = x0.change_ring(A.base_ring())/pnorm(x0, 2)
    Ax = A*x
    vap = float(x * Ax)
    k = 0
    error = tol + 1

    while error > tol and k < itMax:
        xNext = Ax/pnorm(Ax, 2)
        AxNext = A * xNext
        vapNext = float(xNext * AxNext)

        error = abs(vapNext - vap)
        vap = vapNext

        x = xNext
        Ax = AxNext
        k += 1

    return (vap, x, k) 