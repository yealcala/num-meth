# num_meth/methods/poly/sturm.py

"""Provides methods regarding sturm.

- `sign_changes(chain, x)` - Given a sturm chain and a point, counts the number of sign changes in the chain.
- `sturmChain(p)` - Get's the sturm chain given a polynomial.
- `sturmNumRoots(p, evals)` - Counts number of real roots between points of the given list.
- `solveSturm(p, tol)` - Tries to find and approximate all real roots using sturm method.
"""
from sage.all import *
from ...utils import extract_coefficients
from .bound import mclaurin_donut
from ...verify import nonNull

def sign_changes(chain, x: float) -> int:
    """ Given a sturm chain and a point, counts the number of sign changes in the chain.

    Examples: No examples.
        

    Args:
        chain (list): a sturm chain.
        x (float): point where chain will be evaluated.
        

    Returns:
        int: number of sign changes.
    """ 
    nonNull(chain, x)
    values = [RR(p(x)) for p in chain if not p(x).is_zero()]
    
    changes = 0
    for i in range(len(values) - 1):
        if values[i] * values[i+1] < 0:
            changes += 1
    return changes

def sturmChain(p) -> list:
    """ Get's the sturm chain given a polynomial.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(sturmChain(p))
        

    Args:
        p (callable): a polynomial.
        

    Returns:
        list: the sturm chain that corresponds to the polynomial p.
    """ 
    nonNull(p)
    if hasattr(p, 'parent') and hasattr(p.parent(), 'is_poly_ring') and p.parent().is_poly_ring():
        f0 = p
    else:
        coeffs = extract_coefficients(p)
        R = QQ['x']
        f0 = R(list(reversed(coeffs)))

    f1 = f0.derivative()
    chain = [f0, f1]

    while chain[-1].degree() > 0:
        r = -(chain[-2] % chain[-1])
        if r.is_zero():
            break
        chain.append(r)
    return chain

def sturmNumRoots(p, evals: list[float]) -> list[tuple[tuple, int]]:
    """ Counts number of real roots between points of the given list.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(sturmNumRoots(p, [-oo, 0, oo]))
        >>> show(sturmNumRoots(p, [-1.75, -1, 0, 1, 1.75]))

    Args:
        p (callable): a polynomial.
        evals (list): list of numbers.
        

    Returns:
        list: List of tuples containing,
            - 0: a tuple which corresponds to the number of  realroots in an interval (a, b).
            - 1: an integer, the number of real roots in the interval (a, b).
    """ 
    nonNull(p, evals)
    evals = sorted([RR(val) for val in evals])
    if len(evals) < 2:
        raise ValueError("Two points are requiered at least for counting roots in intervals")
    
    chain = sturmChain(p)

    vValues = [sign_changes(chain, pe) for pe in evals]

    out = []
    for i in range(len(evals)-1):
        a, b = evals[i], evals[i+1]
        if a == b:
            continue
            
        nRoots = vValues[i] - vValues[i+1]
        out.append(  ( (a, b), abs(nRoots) )  )

    return out

def solveSturm(p, tol: float = 1e-10) -> tuple[list[float], callable]:
    """ Tries to find and approximate all real roots using sturm method.

    Examples:
        >>> p = 16*x**3 + 12*x**2 - 8*x - 1
        >>> show(solveSturm(p, tol=1e-10))

    Args:
        p (callable): a polynomial.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        

    Returns:
        tuple:
            - 0: list of found real roots.
            - 1: a polynomial which roots are the ones not found (they should be the complex ones, not guaranteed though).
    """ 
    from ..nolinsys.newton import newton
    from ..nolinsys.other import bisection
    nonNull(p, tol)
    _, R_max = mclaurin_donut(p)
    
    chain = sturmChain(p)
    
    def aislarIntervalosReales(a, b):
        n_roots = abs(sign_changes(chain, a) - sign_changes(chain, b))
        if n_roots == 0:
            return []
        if n_roots == 1:
            return [(a, b)]
        
        mid = (a + b) / 2.0
        return aislarIntervalosReales(a, mid) + aislarIntervalosReales(mid, b)

    intervals = aislarIntervalosReales(-R_max, R_max)
    
    raicesReales = []

    xVar = p.variables()[0] if hasattr(p, 'variables') and p.variables() else var('x')
    residuo = p

    for a, b in intervals:
        xBis, _ = bisection(p, (a, b), tol=1e-2)
        
        root, _ = newton(p, xBis, tol=tol)
        raicesReales.append(root)
        
        residuo = (residuo / (xVar - root)).full_simplify()

    return raicesReales, residuo