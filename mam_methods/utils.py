from sage.all import *
from typing import Union

def fEval(f, x_val):
    if hasattr(x_val, '__len__') and not isinstance(x_val, (str, bytes)):
        if hasattr(f, 'variables'):
            vars_f = f.variables()
            subs_dict = {vars_f[i]: x_val[i] for i in range(min(len(vars_f), len(x_val)))}
            return f.subs(subs_dict).change_ring(RR)
        
        if callable(f):
            res = f(x_val)
            return vector(RR, [RR(val).n() for val in res])

    if hasattr(f, 'variables'):
        v_list = f.variables()
        if v_list:
            return RR(f.subs({v_list[0]: x_val})).n()
    
    if callable(f):
        return RR(f(x_val)).n()
    
    return RR(f).n()

def pnorm(v: Union[vector, list, float], p: int) -> float:
    """ Calcula la norma p del vector v con p un entero o infinito (oo) """
    if(p == 0): 
        raise ZeroDivisionError("La norma no puede ser negativa.")
    if type(v) == float:
        return abs(v)
    try:
        iter(v) # type: ignore
    except TypeError:
        return abs(v) # type: ignore
    
    if(p == oo): return max(abs(e) for e in list(v)) # type: ignore
    return sum(abs(e)**p for e in list(v))**(1/p) # type: ignore

def spectral_radius(A: matrix):
    return max([abs(e) for e in list(A.eigenvalues())])

def extract_coefficients(p_input):
    if isinstance(p_input, list):
        return p_input
    
    try:
        if not hasattr(p_input, 'coefficients'):
            R = RR['x']
            p_input = R(p_input)
        
        coeffs = p_input.list()
        coeffs.reverse()
        return coeffs
    except Exception as e:
        raise ValueError(f"No se pudo procesar el polinomio: {e}")