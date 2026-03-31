# num_meth/methods/matrix/other.py

"""Provides some other methods.

The module contains the following functions:

- `companionMatrix(p)` - Gets the companion matrix of a polynomial.
"""
from sage.all import *
from ..poly.all import extract_coefficients

def companionMatrix(p):
    """Gets the companion matrix of a polynomial.

    Examples:
        >>> p = x**4 - 3*x**3 - 60*x**2 + 150*x + 300
        >>> show(companionMatrix(p))
        

    Args:
        p (callable): a polynomial.

    Returns:
        matrix: The companion matrix.
    """
    coeffs = [RR(c) for c in extract_coefficients(p)]
    
    an = coeffs[0]
    coeffs = [c / an for c in coeffs] # Hacemos mónico
    
    n = len(coeffs) - 1
    C = matrix(RR, n, n)
    
    for i in range(1, n):
        C[i, i-1] = 1
        
    for i in range(n):
        C[i, n-1] = -coeffs[n-i]
        
    return C