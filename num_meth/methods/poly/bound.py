# num_meth/methods/poly/bound.py

"""Provides some of the methods for getting and checking bounds of real roots.

- `max_bound_mclaurin(p)` - Calculates the upper mclaurin bound.
- `min_bound_mclaurin(p)` - Calculates the lower mclaurin bound.
- `mclaurin_donut(p)` - Calculates the interval where positive real roots are contained.
- `is_laguerre_upper_positive(p, L)` - Returns whether L is an upper bound of positive real roots.
- `is_laguerre_lower_positive(p, L)` - Returns whether L is a lower bound of positive real roots.
- `is_laguerre_upper_negative(p, L)` - Returns whether L is an upper bound of negative real roots.
- `is_laguerre_lower_negative(p, L)` - Returns whether L is a lower bound of negative real roots.
- `is_newton_upper_positive(p, L)` - Returns whether L is an upper bound of positive real roots.
- `is_newton_lower_positive(p, L)` - Returns whether L is a lower bound of positive real roots.
- `is_newton_upper_negative(p, L)` - Returns whether L is an upper bound of negative real roots.
- `is_newton_lower_negative(p, L)` - Returns whether L is a lower bound of negative real roots.
"""
from sage.all import *
from ...utils import extract_coefficients
from typing import Union
from .eval import horner
from ...verify import nonNull

def max_bound_mclaurin(p):
    """ Calculates the upper mclaurin bound.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(max_bound_mclaurin(p))
        

    Args:
        p (callable): a polynomial.
        

    Returns:
        int: the upper mclaurin bound for real roots.

    Raises:
        ValueError: If first coefficiente is zero.
    """  
    nonNull(p)
    coeffs = extract_coefficients(p)
    n = len(coeffs)

    an = RR(coeffs[0])
    
    if an == 0:
        raise ValueError("First coefficient cannot be zero.")
    
    return max([abs(RR(coeffs[k])/an) for k in range(1, n)])

def min_bound_mclaurin(p):
    """ Calculates the lower mclaurin bound.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(min_bound_mclaurin(p))
        

    Args:
        p (callable): a polynomial.
        

    Returns:
        int: the lower mclaurin bound for real roots.

    Raises:
        ValueError: If first coefficiente is zero.
    """ 
    nonNull(p)
    coeffs = extract_coefficients(p)
    n = len(coeffs)

    a0 = RR(coeffs[n-1])
    
    if a0 == 0:
        raise ValueError("Last coefficient cannot be null.")
    
    return max([abs(RR(coeffs[k])/a0) for k in range(n-1)])

def mclaurin_donut(p):
    """ Calculates the mclaurin interval that contains all real positive roots.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(mclaurin_donut(p))
        

    Args:
        p (callable): a polynomial.
        

    Returns:
        tuple: Represents an interval,
            0: Lower bound.
            1: Upper bound.

    Raises:
        ValueError: If first coefficiente is zero.
    """ 
    nonNull(p)
    return (1/(1 + RR(min_bound_mclaurin(p))), 1 + RR(max_bound_mclaurin(p)))

def is_laguerre_upper_positive(p, L: float) -> bool:
    """ Checks if given number L is an upper bound of real positive roots using Laguerre-Thibaults's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_laguerre_upper_positive(p, 3))
        >>> show(is_laguerre_upper_positive(p, 1.9))
        

    Args:
        p (callable): a polynomial.
        L (float): a positive number.
        

    Returns:
        bool: Whether L is an upper bound.

    Raises:
        ValueError: If L is non-positive.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   

    coeffs = extract_coefficients(p)
    if coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    pVal, Q = horner(coeffs, L)

    return all(c >= 0 for c in Q) and pVal >= 0

def is_laguerre_lower_positive(p, L: float) -> bool:
    """ Checks if given number L is a lower bound of real positive roots using Laguerre-Thibaults's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_laguerre_lower_positive(p, 2.1))
        >>> show(is_laguerre_lower_positive(p, 1))
        

    Args:
        p (callable): a polynomial.
        L (float): a positive number.
        

    Returns:
        bool: Whether L is a lower bound.

    Raises:
        ValueError: If L is non-positive.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   
    
    coeffs = extract_coefficients(p)
    coeffsReverse = list(reversed(coeffs))

    return is_laguerre_upper_positive(coeffsReverse, 1/L)


def reflex_polinomial(coeffs: list) -> list:
    """ P(-x) """

    n = len(coeffs) - 1
    return [coeffs[k] if (n - k) % 2 == 0 else -coeffs[k] for k in range(len(coeffs))]

def is_laguerre_lower_negative(p, L: float) -> bool:
    """ Checks if given number L is a lower bound of real negative roots using Laguerre-Thibaults's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_laguerre_lower_negative(p, -2.1))
        >>> show(is_laguerre_lower_negative(p, -1))
        

    Args:
        p (callable): a polynomial.
        L (float): a negative number.
        

    Returns:
        bool: Whether L is a lower bound.

    Raises:
        ValueError: If L is non-negative.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   
    
    coeffs = extract_coefficients(p)
    coeffsReflejados = reflex_polinomial(coeffs)
    
    return is_laguerre_upper_positive(coeffsReflejados, abs(L))

def is_laguerre_upper_negative(p, L: float) -> bool:
    """ Checks if given number L is an upper bound of real negative roots using Laguerre-Thibaults's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_laguerre_upper_negative(p, -3))
        >>> show(is_laguerre_upper_negative(p, -1.9))
        

    Args:
        p (callable): a polynomial.
        L (float): a negative number.
        

    Returns:
        bool: Whether L is an upper bound.

    Raises:
        ValueError: If L is non-negative.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   
    
    coeffs = extract_coefficients(p)
    coeffsReflejados = reflex_polinomial(coeffs)
    
    coeffsRR = list(reversed(coeffsReflejados))
    
    return is_laguerre_upper_positive(coeffsRR, 1/abs(L))

def is_newton_upper_positive(p, L: float) -> bool:
    """ Checks if given number L is an upper bound of real positive roots using Newton's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_newton_upper_positive(p, 3))
        >>> show(is_newton_upper_positive(p, 1.9))
        

    Args:
        p (callable): a polynomial.
        L (float): a positive number.
        

    Returns:
        bool: Whether L is an upper bound.

    Raises:
        ValueError: If L is non-positive.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   

    coeffs = extract_coefficients(p)
    if coeffs[0] < 0:
        coeffs = [-c for c in coeffs]
        
    n = len(coeffs) - 1
    currentCoeffs = coeffs
    

    for k in range(n + 1):
        pVal, currentCoeffs = horner(currentCoeffs, L)
        
        if pVal < 0:
            return False
        
        if not currentCoeffs and k < n: # Se acaba el grado
            break
            
    return True

def is_newton_lower_positive(p, L: float) -> bool:
    """ Checks if given number L is a lower bound of real positive roots using Newton's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_newton_lower_positive(p, 2.1))
        >>> show(is_newton_lower_positive(p, 1))
        

    Args:
        p (callable): a polynomial.
        L (float): a positive number.
        

    Returns:
        bool: Whether L is a lower bound.

    Raises:
        ValueError: If L is non-positive.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   
    
    coeffs = extract_coefficients(p)    
    coeffsReverse = list(reversed(coeffs))

    return is_newton_upper_positive(coeffsReverse, 1/L)

def is_newton_lower_negative(p, L: float) -> bool:
    """ Checks if given number L is a lower bound of real negative roots using Newton's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_newton_lower_negative(p, -2.1))
        >>> show(is_newton_lower_negative(p, -1))
        

    Args:
        p (callable): a polynomial.
        L (float): a negative number.
        

    Returns:
        bool: Whether L is a lower bound.

    Raises:
        ValueError: If L is non-negative.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   

    coeffs = extract_coefficients(p)
    coeffsReflejados = reflex_polinomial(coeffs)

    return is_newton_upper_positive(coeffsReflejados, abs(L))

def is_newton_upper_negative(p, L: float) -> bool:
    """ Checks if given number L is an upper bound of real negative roots using Newton's method.

    Examples:
        >>> p = x**2 - 4
        >>> show(is_newton_upper_negative(p, -3))
        >>> show(is_newton_upper_negative(p, -1.9))
        

    Args:
        p (callable): a polynomial.
        L (float): a negative number.
        

    Returns:
        bool: Whether L is an upper bound.

    Raises:
        ValueError: If L is non-negative.
    """ 
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   
    
    coeffs = extract_coefficients(p)
    coeffsReflejados = reflex_polinomial(coeffs)
    coeffsRR = list(reversed(coeffsReflejados))

    return is_newton_upper_positive(coeffsRR, 1/abs(L))