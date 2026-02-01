from sage.all import *
from ...utils import extract_coefficients
from ..nolinsys.newton import newtonSystem
from ..poly.eval import horner
from .bound import reflejar_polinomio
from ...verify import nonNull

def __sign_changes(lista) -> int:
    filtered = [x for x in lista if x != 0]
    if len(filtered) < 2: return 0

    return sum(1 for i in range(len(filtered) - 1) 
                 if filtered[i] * filtered[i+1] < 0)

def descartes(p: callable) -> tuple[int, int]:
    nonNull(p)
    coeffs = [RR(c) for c in extract_coefficients(p)]
    coeffs_reflejados = reflejar_polinomio(coeffs)
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
def solveBairstow(p, u0: vector, v0: vector, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[any, any, tuple[float, float]]:
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