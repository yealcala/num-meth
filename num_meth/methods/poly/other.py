# num_meth/methods/poly/other.py

"""Provides some other methods.

- `descartes(p)` - Return the possible max number of real positive and negative roots.
- `solveBairstow(p, u0, v0, tol, itMax, minIt)` - Applies bairstow method for deflacting polynomial using complex roots.
- `graeffe(p, it)` - Finds all real roots using Graeffe method.

"""
from sage.all import *
from ...utils import extract_coefficients
from ..nolinsys.newton import newtonSystem
from ..poly.eval import horner
from .bound import reflex_polinomial
from ...verify import nonNull

def __sign_changes(lista) -> int:
    filtered = [x for x in lista if x != 0]
    if len(filtered) < 2: return 0

    return sum(1 for i in range(len(filtered) - 1) 
                 if filtered[i] * filtered[i+1] < 0)

def descartes(p: callable) -> tuple[int, int]:
    """ Get's maximum possible real positive and negative roots by counting sign changes in coefficients.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(descartes(p))
        

    Args:
        p (callable): a polynomial.
        

    Returns:
        tuple:
            - 0: maximum number of real positive roots.
            - 1: maximum number of real negative roots.

    Raises:
        ValueError: If polynomial cannot be processed.
    """ 
    nonNull(p)
    coeffs = [RR(c) for c in extract_coefficients(p)]
    coeffs_reflejados = reflex_polinomial(coeffs)
    return (__sign_changes(coeffs), __sign_changes(coeffs_reflejados))


def __bairstow_residuals(params, coeffs):
    u, v = params[0], params[1]
    n = len(coeffs) - 1
    
    b = [0.] * (n + 1)

    b[0] = coeffs[0]
    b[1] = coeffs[1] + u * b[0]
    for i in range(2, n + 1):
        b[i] = coeffs[i] + u * b[i-1] + v * b[i-2]


    
    return vector([b[n], b[n-1]])
def solveBairstow(p, u0: float, v0: float, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[any, any, tuple[float, float]]:
    """ Searches for complex roots applying bairstow method for some initial x^2 - u0*x - v0

    Examples:
        >>> p = x**3 - x - 1
        >>> solveBairstow(p, 1.1, 0.9, 1e-10)
        

    Args:
        p (callable): a polynomial.
        u0 (float): starting number of x^2 - u0*x - v0
        v0 (float): starting number of x^2 - u0*x - v0
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        minIt (Optional[int]): Minimum iterations. Default: 0.
        

    Returns:
        tuple:
            - 0: a complex root.
            - 1: the other complex root (conjugated)
            - 2: tuple:
                - 0: Approximated coefficient u
                - 1: Approximated coefficient v

    Raises:
        ValueError: If polynomial cannot be processed.
        ConvergenceError: If divergence has been found.
    """ 
    nonNull(p, u0, v0, tol, itMax, minIt)
    coeffs = [RR(c) for c in extract_coefficients(p)]

    F = lambda params: __bairstow_residuals(params, coeffs)
    x0 = vector(RR, [u0, v0])
    res, _ = newtonSystem(F, x0, tol=tol, itMax=itMax, minIt=minIt)
    u = RR(res[0])
    v  = RR(res[1])
    a = u/2
    b = sqrt(abs(-v - u**2/4))
    
    return a+b*I, a-b*I, (u, v)

def __graeffeSeparate(coeffs: list) -> list:
    n = len(coeffs) - 1
    new_coeffs = [0.0] * (n + 1)
    
    for k in range(n + 1):
        sum_terms = coeffs[k]**2
        
        j = 1
        while (k - j >= 0) and (k + j <= n):
            term = 2 * (-1)**j * coeffs[k-j] * coeffs[k+j]
            sum_terms += term
            j += 1
            
        new_coeffs[k] = sum_terms
        
    return new_coeffs

def graeffe(p, it: int = 5) -> list:
    """ Searches roots by separating and approximating them.

    Examples:
        >>> p = x**3 - 4*x**2 - 7*x + 10
        >>> roots = graeffe(p, 5)
        

    Args:
        p (callable): a polynomial.
        it (int): number of iterations (= times roots will be separated)
        

    Returns:
        list: list of all roots found this way.

    Raises:
        ValueError: If polynomial cannot be processed.
    """ 
    nonNull(p, it)
    coeffs = [RR(c) for c in extract_coefficients(p)]
    n = len(coeffs) - 1
    
    current_coeffs = coeffs
    for _ in range(it):
        current_coeffs = __graeffeSeparate(current_coeffs)
    
    N = 2**it
    roots_approx = []
    

    for i in range(1, n + 1):
        try:
            roots_approx.append(abs(current_coeffs[i] / current_coeffs[i-1])**(1/N))
        except ZeroDivisionError:
            roots_approx.append(0.0)


    final_roots = []
    for approx in roots_approx: # Poner signos
        if approx < 1e-15:
            final_roots.append(0.0)
            continue
            
        p_pos, _ = horner(coeffs, approx)
        p_neg, _ = horner(coeffs, -approx)
        
        if abs(p_pos) <= abs(p_neg):
            final_roots.append(approx)
        else:
            final_roots.append(-approx)
            
    return sorted(final_roots)