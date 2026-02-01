from sage.all import *
from ...verify import matrixMustBeSquare

def leverrier(A: matrix):
    """
    Calcula los coeficientes del polinomio característico p(x) = x^n + c1*x^(n-1) + ... + cn
    basado en las trazas de las potencias de la matriz.
    """
    matrixMustBeSquare(A)
    n = A.ncols()
    base_ring = A.base_ring()
    R = PolynomialRing(base_ring, 'x')

    s = []
    p = []
    
    Ak = A
    
    for k in range(1, n + 1):
        sk = Ak.trace()
        s.append(sk)
        
        suma = sk
        for j in range(1, k):
            suma += p[j-1] * s[k-j-1]
        
        pk = -1/k * suma
        p.append(pk)
        
        if k < n:
            Ak = Ak * A

    return R(p[::-1] + [base_ring(1)] )
