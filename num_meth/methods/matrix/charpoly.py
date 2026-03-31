# num_meth/methods/matrix/charpoly.py

"""Provides methods for characteristic polynomials.

The module contains the following functions:

- `leverrier(A)` - Returns the characteristic polynomial of a matrix using Leverrier method.
"""
from sage.all import *
from ...verify import matrixMustBeSquare

def leverrier(A: matrix):
    """Calculates the characteristic polynomial of a matrix using Leverrier method.

    Examples:
        >>> A = matrix(QQ, [[3, 1, 5], [3, 3, 1], [4, 6, 4]])
        >>> show("Leverrier:", leverrier(A))
        

    Args:
        A (matrix): A square matrix

    Returns:
        callable: Returns the characteristic polynomial of the matrix.
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
