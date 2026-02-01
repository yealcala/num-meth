from sage.all import *
from ..matrix.power import powerIteration
from ...utils import pnorm
from typing import Union
from ...verify import nonNull, matrixMustBeSquare

def wielandt(A: matrix, vap: float, vep: vector, method: Union[str, None] = "first") -> tuple[matrix, matrix, vector, int]:
    nonNull(A, vap, vep)
    matrixMustBeSquare(A)
    n = A.ncols()
    if A.base_ring() == QQ or A.base_ring() == ZZ:
        A = A.change_ring(RR)

    i = 0
    method = method.lower() if method is not None else method

    if method == "partial":
            maxVal = -1
            for imax in range(n):
                if abs(vep[imax]) > maxVal:
                    maxVal = abs(vep[imax])
                    i = imax
    else:
        for imin in range(n):
            if vep[imin] != 0: 
                i = imin
                break
        
    X = 1/(vap*vep[i])*A.row(i)
    B = A - vap*vep.column() * X.row()
    Bred = B.delete_rows([i]).delete_columns([i])
    return (B, Bred, X, i)

def wielandtSpectrum(A: matrix, tol: float = 1e-10, itMax: int = 100, method: Union[str, None] = "first") -> list:
    nonNull(A, tol, itMax)
    n = A.ncols()
    if n == 1: return [(A[0,0], vector(A.base_ring(), [1]))]
    vap1, vep1, _ = powerIteration(A, vector(A.base_ring(), [1]*n), tol, itMax)
    out = [(vap1, vep1)]
    _, Bred, X, imin = wielandt(A, vap1, vep1, method)
    listaRed = wielandtSpectrum(Bred, tol, itMax)
    for tupla in listaRed:
        vap, vepRed = tupla[0], tupla[1]
        wiList = list(vepRed)
        wiList.insert(imin, 0)
        wi = vector(A.base_ring(), wiList)
        vep = (vap - vap1)*wi + vap1* (X*wi) *vep1 
        out.append((vap, vep/pnorm(vep, oo)))
    return out

def householder(A: matrix, vep: vector) -> tuple[matrix, matrix, matrix, float]:
    nonNull(A, vep)
    matrixMustBeSquare(A)

    n = A.ncols()
    base_ring = A.base_ring()
    
    s = pnorm(vep, 2) * (-1 if vep[0] > 0 else 1)
    
    e1 = vector(base_ring, [1] + [0]*(n-1))
    u = vep - s * e1
    
    uu = u * u
    if uu == 0:
        P = matrix.identity(base_ring, n)
    else:
        P = matrix.identity(base_ring, n) - (2/uu) * u.column() * u.row()
    

    Ap = P * A * P
    
    ct = Ap.submatrix(0, 1, 1, n-1)
    Atilde = Ap.delete_rows([0]).delete_columns([0])
    
    return (Ap, Atilde, P, ct)

def householderSpectrum(A: matrix, tol: float = 1e-10, itMax: int = 100) -> list:
    nonNull(A, tol, itMax)
    matrixMustBeSquare(A)
    n = A.ncols()
    
    if A.base_ring() == QQ or A.base_ring() == ZZ:
        A = A.change_ring(RR)
    
    base_ring = A.base_ring()
    
    if n == 1: 
        return [(A[0,0], vector(base_ring, [1]))]
    
    vap1, vep1, _ = powerIteration(A, vector(base_ring, [1]*n), tol, itMax)
    out = [(vap1, vep1)]
    
    _, Atilde, P, ct = householder(A, vep1)
    
    listaRed = householderSpectrum(Atilde, tol, itMax)
    
    for mu, v_tilde in listaRed:
        if abs(mu - vap1) < tol:
            b = 0
        else:
            b = (ct * v_tilde.column())[0,0] / (mu - vap1)
        
        u_prime = vector(base_ring, [b] + list(v_tilde))
        
        vep = P * u_prime
        
        out.append((mu, vep / pnorm(vep, oo)))
        
    return out