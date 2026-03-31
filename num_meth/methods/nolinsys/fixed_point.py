# num_meth/methods/nolinsys/fixed_point.py

"""Provides some of the methods for fixed point

- `fixedPointIteration(g, x0, tol, itMax)` - Applies fixed point iteration.
- `fixedPointBanach(g, interval, tol, itMax, minIt)` - Applies fixed point iteration checking Banach theorem.
- `fixedPointSystem(G, x0, tol, itMax, minIt)` - Applies fixed point iteration for a given system.
- `schroeder(g, p, maxP, tol)` - Convergence order using schroeder.
"""

from sage.all import *
from sage.modules.free_module_element import FreeModuleElement
from ...utils import pnorm, fEval
from ...exceptions import ConvergenceError
from ..matrix.power import powerIteration
from .other import findGeneralRoots
from typing import Union
from ...verify import nonNull

def fixedPointIteration(g: callable, x0: Union[vector, float], tol: float = 1e-10, 
                        itMax: int = 100, minIt: int = 0) -> tuple[Union[float, vector], int]:
    """ Iterates over x_{n+1} = g(x_n)

    Examples:
        >>> g = arctan(4-x)
        >>> show(fixedPointIteration(g, 1.5, 1e-5, minIt=3))
        

    Args:
        g (callable): any callable function.
        x0 (float): Starting point of iterations.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        minIt (Optional[int]): Minimum iterations. Default: 0.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
    """  
    nonNull(g, x0, tol, itMax, minIt)
    xCurr = RR(x0) if not isinstance(x0, (tuple, list, FreeModuleElement)) else x0
    for k in range(1, itMax + 1):
        xNext = fEval(g, xCurr)
        
        if pnorm(xNext - xCurr, 2) < tol and k >= minIt:
            return xNext, k
            
        xCurr = xNext
    return xCurr, k

def fixedPointBanach(g: callable, interval: tuple[float, float], 
                     tol: float = 1e-10, itMax: int = 100, minIt: int = 0):
    """ Iterates over x_{n+1} = g(x_n) checking first if converges by Banach.

    Examples:
        >>> g = arctan(4-x)
        >>> show(fixedPointBanach(g, (1, 2), 1e-5, minIt=3))
        

    Args:
        g (callable): any callable function.
        interval (tuple[float, float]): Interval of floats.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        minIt (Optional[int]): Minimum iterations. Default: 0.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
    """  
    nonNull(g, interval, tol, itMax, minIt)
    a, b = interval
    if a >= b:
        raise ValueError("Interval is not valid!")
    x_var = g.variables()[0]
    
    dg = diff(g, x_var)
    puntos_criticos_g = findGeneralRoots(dg, interval)

    # --- ESTUDIO DEL INTERVALO ---
    candidatos_extremos = [a, b] + [r for r in puntos_criticos_g if a < r < b]
    valores_g = [g.subs({x_var: c}).n() for c in candidatos_extremos]
    
    g_min, g_max = min(valores_g), max(valores_g)
    
    if g_min < a or g_max > b:
        raise ValueError(
            f"Fallo de Auto-mapeo: g([{a}, {b}]) = [{g_min:.4f}, {g_max:.4f}]. "
            f"No se garantiza que el punto fijo esté en el intervalo."
        )

    # --- ESTUDIO DE LA CONTRACCIÓN ---
    ddg = diff(dg, x_var)
    if ddg.is_constant():
        puntosCriticosDg = []
    else:
        puntosCriticosDg = findGeneralRoots(ddg, interval)
        
    candidatos_L = [a, b] + [r for r in puntosCriticosDg if a < r < b]
    l_values = [abs(dg.subs({x_var: c}).n()) for c in candidatos_L]
    L = max(l_values)
    
    if L >= 1:
        raise ConvergenceError(f"La función no es contractiva en el intervalo. L = {L:.4f} >= 1.")

    return fixedPointIteration(g, (interval[0] + interval[1]) / 2, 
                               tol, itMax, minIt)


def fixedPointSystem(G: callable, x0: vector, tol: float = 1e-10, itMax: int = 100, minIt: int = 0):
    """ Iterates over a system x_{n+1} = g(x_n)

    Examples:
        >>> def G(v):
        >>>     return vector([ sin(v[0] + v[1]) - v[0], cos(v[0] - v[1]) - v[1] ])

        >>> x0 = vector(RR, [1, 1])
        >>> sol, k = fixedPointSystem(G, x0, tol=1e-10)
        >>> show(f"Solution: {sol}")
        

    Args:
        G (callable): a function as seen in the example.
        x0 (vector): Starting point of iterations.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 100.
        minIt (Optional[int]): Minimum iterations. Default: 0.
        

    Returns:
        tuple: A tuple with two components,
            - 0: solution after iterations.
            - 1: Number of iterations made.

    Raises:
        ValueError: If some argument is not valid.
        ConvergenceError: If it diverges at least in first iteration.
    """  
    nonNull(G, x0, tol, itMax, minIt)
    n = len(x0)
    

    vars = [var(f'v{i}') for i in range(n)]
    G_sym = G(vector(vars))
    
    J = jacobian(G_sym, vars)
    

    subs_dict = {vars[i]: x0[i] for i in range(n)}
    J0 = J.subs(subs_dict).change_ring(RR)
    
    rho, _, _ = powerIteration(J0, vector(RR, [1]*n), tol=1e-5)
    
    if rho >= 1:
        raise ConvergenceError(f"No se garantiza convergencia: rho(J_G(x0)) = {rho:.4f} >= 1")
    
    # Si pasa el control, procedemos al punto fijo escalar/vectorial
    return fixedPointIteration(G, x0, tol, itMax, minIt)

def schroeder(g: callable, p: float, itMax: int = 10, tol: float = 1e-12) -> tuple[Union[int, None], Union[float, None]]:
    """ Calculates convergence order using Schroeder's theorem.

    Examples:
        >>> g = lambda x: (x**2 + 2) / (2*x)
        >>> fixedPoint = sqrt(2).n()

        >>> show(schroeder(g, fixedPoint))
        

    Args:
        g (callable): any callable function.
        p (vector): A fixed point.
        tol (Optional[float]): Maximum allowed tolerance. Default: 10^{-10}.
        iMax (Optional[int]): Maximum iterations allowed. Default: 10.
        

    Returns:
        tuple: A tuple with two components,
            - 0: Order of convergence or None if itMax was reached.
            - 1: Asympthotic value or None if itMax was reached.

    Raises:
        ValueError: If some argument is not valid.
    """
    nonNull(g, p, itMax, tol)
    t = var('t')
    gSym = g(t)

    if abs(p - g(p)) >= tol:
        raise ValueError(f"El valor {p} no es un punto fijo válido de g bajo la tolerancia {tol}")
    
    for k in range(1, itMax):
        dgSym = diff(gSym, t, k)
        
        deriv = dgSym.subs(t=p).n()
        
        if abs(deriv) > tol:
            return (k, abs(deriv)/factorial(k))
            
    return (None, None)