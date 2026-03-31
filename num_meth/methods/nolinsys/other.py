# num_meth/methods/nolinsys/other.py

"""Provides other methods related to non linear systems.

- `bisection(f, interval, tol, itMax)` - Applies bisection iteration.
- `regulaFalsi(f, interval, tol, itMax)` - Applies Regula Falsi iteration.
- `secantIteration(f, x0, x1, tol, itMax)` - Applies Secant iteration.
- `mullerIteration(P, x0, x1, x2, tol, itMax, minIt)` - Applies Müller iteration.
- `mullerIteration(P, x0, x1, x2, tol, itMax, minIt)` - Applies Müller iteration.
- `findGeneralRoots(f, interval, steps, tol)` - Tries to find all real roots.
"""

from sage.all import *

from ...exceptions import ConvergenceError
from ...utils import fEval
from ...verify import nonNull


def bisection(f: callable, interval: tuple[float, float], tol: float = 1e-10, itMax: int = 100) -> tuple[float, int]:
    """ Bisection method. Finds root by reducing interval's length.

    Examples:
        >>> f = sqrt(x)*sin(x) - x**3 + 2
        >>> show(bisection(f, (1, 2), tol=1/30))
        

    Args:
        f (callable): any callable function.
        interval (tuple[float, float]): an interval.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
        ConvergenceError: If convergence is not guaranteed.
    """  
    nonNull(f, interval, tol, itMax)
    a, b = interval
    if a >= b:
        raise ValueError("Interval is not valid!")
    fa = fEval(f, a)
    fb = fEval(f, b)
    
    if fa * fb > 0:
        raise ConvergenceError(f"Interval ({a}, {b}) doesn't guarantee any root (f(a)*f(b) > 0)!")
    
    for k in range(1, itMax + 1):
        c = a + (b-a)/2 # Más estable
        fc = fEval(f, c)

        if (b-a)/2 < tol or fc == 0:
            return c, k
        
        if fa*fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return (a+b)/2, itMax

def regulaFalsi(f: callable, interval: tuple[float, float], tol: float = 1e-10, itMax: int = 100):
    """ Regula falsi method.

    Examples:
        >>> f = sqrt(x)*sin(x) - x**3 + 2
        >>> show(regulaFalsi(f, (1, 2), tol=1/30))
        

    Args:
        f (callable): any callable function.
        interval (tuple[float, float]): an interval.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
        ConvergenceError: If convergence is not guaranteed.
    """  
    nonNull(f, interval, tol, itMax)
    a, b = RR(interval[0]), RR(interval[1])
    fa, fb = fEval(f, a), fEval(f, b)
    if a >= b:
        raise ValueError("El intervalo es inválido!")
    elif fa * fb > 0:
        raise ConvergenceError(f"Interval ({a}, {b}) doesn't guarantee any root (f(a)*f(b) > 0)!")
    
    cPrev = a

    for k in range(1, itMax + 1):
        c = b - fb*(b-a)/(fb - fa)
        fc = fEval(f, c)

        if abs(c - cPrev) < tol or abs(fc) < tol:
            return c, k
        cPrev = c

        if fa*fc < 0:
            b = c
            fb = fc
        else:
            a = c
            fa = fc
    return c, itMax


def secantIteration(f: callable, x0: float, x1: float, tol: float = 1e-10, itMax: int = 100) -> tuple[float, int]:
    """ Secant iteration method.

    Examples:
        >>> ec = x**6 - x - 1

        >>> show(secantIteration(ec, 1, 1.5, 5e-5))
        

    Args:
        f (callable): any callable function.
        x0 (float): first point of iteration.
        x1 (float): second point of iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
        ConvergenceError: If convergence is not guaranteed.
    """  
    nonNull(f, x0, x1, tol, itMax)
    f0 = fEval(f, x0)
    f1 = fEval(f, x1)
    
    if f0 * f1 > 0:
        raise ConvergenceError("Can't guarantee any root for these initial points!")

    for k in range(1, itMax + 1):
        if abs(f1 - f0) < tol:
            raise ZeroDivisionError("Difference of images is null!.")
            
        xNext = x1 - f1 * (x1 - x0) / (f1 - f0)
        fNext = fEval(f, xNext)
        
        ek = abs(xNext - x1)
        
        if ek < tol:
            return xNext, k
            
        x0, f0 = x1, f1
        x1, f1 = xNext, fNext
        
    return x1, itMax

def mullerIteration(P: callable, x0: float, x1: float, x2: float, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[float, int]:
    """ Muller iteration method. Approximates root using parabolas.

    Examples:
        >>> f = lambda x: cos(x) - x**3
        >>> show("Muller:", mullerIteration(f, 0, 1, 0.5, tol=1e-10, itMax=100))
        

    Args:
        f (callable): any callable function.
        x0 (float): first point of iteration.
        x1 (float): second point of iteration.
        x2 (float): third point of iteration.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        minIt (Optional[int]): Minimum of iterations. Default: 0.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
        ConvergenceError: If divergence is found.
    """  
    nonNull(P, x0, x1, x2, tol, itMax, minIt)
    x0, x1, x2 = RR(x0), RR(x1), RR(x2)
    
    f = lambda val: fEval(P, val)

    for i in range(itMax):
        # Evaluaciones en los tres puntos actuales
        f0, f1, f2 = f(x0), f(x1), f(x2)
        
        d02 = x0 - x2
        d12 = x1 - x2
        d01 = x0 - x1
        f12 = f1 - f2
        

        denom = d02 * d12 * d01
        a = (d12 * (f0 - f2) - d02 * f12) / denom
        b = (d02**2 * f12 - d12**2 * (f0 - f2)) / denom
        c = f2
        

        discriminante = sqrt(b**2 - 4*a*c)
        if abs(b + discriminante) > abs(b - discriminante):
            den = b + discriminante
        else:
            den = b - discriminante
            
        dx = -2 * c / den
        x3 = x2 + dx
        
        # Criterio de parada con periodo de gracia minIt
        if i >= minIt and abs(dx) < tol:
            return x3, i
            
        # Desplazamiento de puntos para la siguiente iteración
        x0, x1, x2 = x1, x2, x3
        
    raise ConvergenceError(f"Müller couldn't converge after {itMax} iterations.")

def findGeneralRoots(f: callable, interval: tuple[float, float], steps: int = 200, tol: float = 1e-10) -> list:
    """ This method is not intended for public use, so use it at your own risk!
     
        Tried to approximate general real roots by dividing interval in subintervals and applying both bisection and Newton to each interval if possible.

    Examples: No examples.
        

    Args:
        f (callable): any callable function.
        interval (tuple[float, float]): an interval.
        steps (Optional[int]): steps. Default: 200.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        

    Returns:
        tuple: A list with possible real roots. May be all real roots in the interval, may be none.

    Raises:
        ValueError: If some argument is not valid.
    """
    from .newton import newton
    nonNull(f, interval, steps, tol)

    a, b = interval
    if a >= b:
        raise ValueError("Interval is not valid!")

    xVar = f.variables()[0]
    dfSym = diff(f, xVar)
    
    HSym = f / dfSym
    
    fNum = lambda val: RR(f.subs({xVar: val}).n())
    HNum = lambda val: RR(HSym.subs({xVar: val}).n())
    
    h = (b - a) / steps
    roots = []
    
    for i in range(steps):
        izq = a + i*h
        dcha = izq + h
        
        try:
            h_l, h_r = HNum(izq), HNum(dcha)
            
            if h_l * h_r <= 0:
                rootApprox, _ = bisection(HNum, (izq, dcha), tol=1e-2)
                rootFinal, _ = newton(fNum, rootApprox, tol=tol)
                roots.append(rootFinal)
        except ZeroDivisionError:
            roots.append(RR(izq))
            
    return sorted(list(set(roots)))