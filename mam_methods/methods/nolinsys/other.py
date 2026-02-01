from sage.all import *

from ...exceptions import ConvergenceError
from ...utils import fEval
from ...verify import nonNull


def bisection(f: callable, interval: tuple[float, float], tol: float = 1e-10, itMax: int = 100) -> tuple[float, int]:
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
    """
    Método de la Secante. Requiere dos puntos iniciales y no utiliza derivadas.
    Mantiene la consistencia con la arquitectura de la librería.
    """
    nonNull(f, x0, x1, tol, itMax)
    f0 = fEval(f, x0)
    f1 = fEval(f, x1)
    
    if f0 * f1 > 0:
        raise ConvergenceError("Can't guarantee any root for these initial points!")

    for k in range(1, itMax + 1):
        if abs(f1 - f0) < tol:
            raise ZeroDivisionError("Diferencia de funciones nula: Secante fallida.")
            
        xNext = x1 - f1 * (x1 - x0) / (f1 - f0)
        fNext = fEval(f, xNext)
        
        ek = abs(xNext - x1)
        
        if ek < tol:
            return xNext, k
            
        x0, f0 = x1, f1
        x1, f1 = xNext, fNext
        
    return x1, itMax

def mullerIteration(P: callable, x0: float, x1: float, x2: float, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[float, int]:
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

def findGeneralRoots(f: callable, interval: tuple[int, int], steps: int = 200, tol: float = 1e-10) -> list:
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