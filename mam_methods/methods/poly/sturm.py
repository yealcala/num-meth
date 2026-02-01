from sage.all import *
from ...utils import extract_coefficients
from .bound import mclaurin_donut
from ...verify import nonNull

def sign_changes(chain, x: float) -> int:
    nonNull(chain, x)
    values = [RR(p(x)) for p in chain if not p(x).is_zero()]
    
    changes = 0
    for i in range(len(values) - 1):
        if values[i] * values[i+1] < 0:
            changes += 1
    return changes

def sturmChain(p) -> list:
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
    nonNull(p, evals)
    evals = sorted([RR(val) for val in evals])
    if len(evals) < 2:
        raise ValueError("Se requieren al menos dos puntos en la lista.")
    
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