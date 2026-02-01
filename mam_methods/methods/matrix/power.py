from sage.all import *
from ...utils import pnorm
from ..linsys.factorization import doolittlePLU
from ..linsys.remonte import remonteArriba, remonteAbajo
from ...verify import nonNull, matrixMustBeSquare, matrixMustBeSymmetric, vectorMustBeDimension


def powerIteration(A: matrix, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo):
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

def shiftedPowerIteration(A: matrix, d: float, z0, tol: float = 1e-10, itMax: int = 100):
    nonNull(A, d, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())
    n = A.nrows()

    As = A - d * matrix.identity(A.base_ring(), n)
    mu, vep, k = powerIteration(As, z0, tol, itMax)
    
    return (mu + d, vep, k)

def inversePowerIteration(A: matrix, z0, tol: float = 1e-10, itMax: int = 100, norm: int = oo):
    nonNull(A, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())
    z = z0.change_ring(A.base_ring()) / pnorm(z0, norm)

    P, L, U = doolittlePLU(A) 
    
    lOld = 0
    k = 0
    error = tol + 1
    while error > tol and k < itMax:
        y = remonteArriba(U, remonteAbajo(L, P*z))
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

def shiftedInversePowerIteration(A: matrix, d: float, z0, tol: float = 1e-10, itMax: int = 100):
    nonNull(A, d, z0, itMax, norm)
    matrixMustBeSquare(A)
    vectorMustBeDimension(z0, A.ncols())
    n = A.nrows()

    As = A - d * matrix.identity(A.base_ring(), n)
    mu, vep, k = inversePowerIteration(As, z0, tol, itMax)
    
    return (mu + d, vep, k)

def rayleigh(A: matrix, x0, tol: float = 1e-10, itMax: int = 100):
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