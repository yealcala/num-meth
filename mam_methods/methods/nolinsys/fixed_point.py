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

def schroeder(g: callable, p: float, maxP: int = 10, tol: float = 1e-12) -> tuple[Union[int, None], Union[float, None]]:
    nonNull(g, p, maxP, tol)
    t = var('t')
    gSym = g(t)

    if abs(p - g(p)) >= tol:
        raise ValueError(f"El valor {p} no es un punto fijo válido de g bajo la tolerancia {tol}")
    
    for k in range(1, maxP):
        dgSym = diff(gSym, t, k)
        
        deriv = dgSym.subs(t=p).n()
        
        if abs(deriv) > tol:
            return (k, abs(deriv)/factorial(k))
            
    return (None, None)