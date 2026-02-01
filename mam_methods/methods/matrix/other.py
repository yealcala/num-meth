from sage.all import *
from ..poly.all import extract_coefficients

def companionMatrix(p):
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