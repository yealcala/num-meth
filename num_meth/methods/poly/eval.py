# num_meth/methods/poly/eval.py

"""Provides some of the methods for evaluating polynomials.

- `horner(p, x)` - Evaluates polynomial on a point.
"""

from sage.all import *
from ...utils import extract_coefficients
from ...verify import nonNull

def horner(p, x: float) -> tuple[float, callable]:
    """ Evaluates polynomial using Horner (Ruffini) method.

    Examples:
        >>> p = x**2 - 4
        >>> show(horner(p, 10))
        

    Args:
        p (callable): a polynomial.
        x (float): a number.
        

    Returns:
        tuple:
            - 0: a number which corresponds to evaluating polynomial, p(x)
            - 1: list of coefficients that results when making deflaction with (r-x).

    Raises:
        ValueError: If first coefficiente is zero.
    """ 
    nonNull(p, x)
    coeffs = extract_coefficients(p)
    n = len(coeffs) 
    x = RR(x)
    

    Q = [RR(coeffs[0])]

    pVal = RR(coeffs[0])
    
    for i in range(1, n):
            
        pVal = RR(coeffs[i]) + x * pVal
        if i < n - 1:
            Q.append(pVal)


    if n == 1:
        return pVal, []
        
    return pVal, Q