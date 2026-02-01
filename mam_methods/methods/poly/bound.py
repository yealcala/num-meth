from sage.all import *
from ...utils import extract_coefficients
from typing import Union
from .eval import horner
from ...verify import nonNull

def max_bound_mclaurin(p):
    nonNull(p)
    coeffs = extract_coefficients(p)
    n = len(coeffs)

    an = RR(coeffs[0])
    
    if an == 0:
        raise ValueError("First coefficient cannot be null.")
    
    return max([abs(RR(coeffs[k])/an) for k in range(1, n)])

def min_bound_mclaurin(p):
    nonNull(p)
    coeffs = extract_coefficients(p)
    n = len(coeffs)

    a0 = RR(coeffs[n-1])
    
    if a0 == 0:
        raise ValueError("Last coefficient cannot be null.")
    
    return max([abs(RR(coeffs[k])/a0) for k in range(n-1)])

def mclaurin_donut(p):
    nonNull(p)
    return (1/(1 + RR(min_bound_mclaurin(p))), 1 + RR(max_bound_mclaurin(p)))

def is_laguerre_superior_positive(p, L: float) -> Union[bool, ValueError]:
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   

    coeffs = extract_coefficients(p)
    if coeffs[0] < 0:
        coeffs = [-c for c in coeffs]

    pVal, Q = horner(coeffs, L)

    return all(c >= 0 for c in Q) and pVal >= 0

def is_laguerre_inferior_positive(p, L: float) -> Union[bool, ValueError]:
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   
    
    coeffs = extract_coefficients(p)
    coeffsReverse = list(reversed(coeffs))

    return is_laguerre_superior_positive(coeffsReverse, 1/L)


def reflejar_polinomio(coeffs: list) -> list:
    """ P(-x) """

    n = len(coeffs) - 1
    return [coeffs[k] if (n - k) % 2 == 0 else -coeffs[k] for k in range(len(coeffs))]

def is_laguerre_inferior_negative(p, L: float) -> Union[bool, ValueError]:
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   
    
    coeffs = extract_coefficients(p)
    coeffsReflejados = reflejar_polinomio(coeffs)
    
    return is_laguerre_superior_positive(coeffsReflejados, abs(L))

def is_laguerre_superior_negative(p, L: float) -> Union[bool, ValueError]:
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   
    
    coeffs = extract_coefficients(p)
    coeffsReflejados = reflejar_polinomio(coeffs)
    
    coeffsRR = list(reversed(coeffsReflejados))
    
    return is_laguerre_superior_positive(coeffsRR, 1/abs(L))

def is_newton_superior_positive(p, L: float) -> bool:
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

def is_newton_inferior_positive(p, L: float) -> bool:
    """ Cota inferior positiva usando el polinomio recíproco y Newton. """
    nonNull(p, L)
    L = RR(L)
    if L <= 0:
        raise ValueError("Bound to be checked must be positive!")   
    
    coeffs = extract_coefficients(p)    
    coeffsReverse = list(reversed(coeffs))

    return is_newton_superior_positive(coeffsReverse, 1/L)

def is_newton_inferior_negative(p, L: float) -> bool:
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   

    coeffs = extract_coefficients(p)
    coeffsReflejados = reflejar_polinomio(coeffs)

    return is_newton_superior_positive(coeffsReflejados, abs(L))

def is_newton_superior_negative(p, L: float) -> bool:
    nonNull(p, L)
    L = RR(L)
    if L >= 0:
        raise ValueError("Bound to be checked must be negative!")   
    
    coeffs = extract_coefficients(p)
    coeffsReflejados = reflejar_polinomio(coeffs)
    coeffsRR = list(reversed(coeffsReflejados))

    return is_newton_superior_positive(coeffsRR, 1/abs(L))