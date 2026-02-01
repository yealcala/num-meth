from sage.all import *
from .fixed_point import *
from ..linsys.all import solvePLU
from ...verify import nonNull

def __fEval(f, x_val):
    if hasattr(f, 'variables'):
        v_list = f.variables()
        if v_list:
            return f.subs({v_list[0]: x_val})
    
    if callable(f):
        return f(x_val)
    
    try:
        return f
    except:
        return f(x_val)

def newtonFourierSemilocal(f: callable, interval: tuple[float, float], x0: float, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[float, float]:
    from .other import findGeneralRoots
    nonNull(f, interval, x0, tol, itMax, minIt)
    a, b = interval
    if a >= b:
        raise ValueError("Interval is not valid!")
    x_var = f.variables()[0]
    
    dfSym = diff(f, x_var)
    ddf_sym = diff(dfSym, x_var)
    

    fa = f.subs({x_var: a}).n()
    fb = f.subs({x_var: b}).n()
    if fa * fb >= 0:
        raise ConvergenceError("Bolzano no se cumple: f(a) y f(b) deben tener signos opuestos.")


    raices_df = findGeneralRoots(dfSym, interval)
    if any(a < r < b for r in raices_df):
        raise ConvergenceError("f'(x) no puede tener raíces en el intervalo.")



    raices_ddf = findGeneralRoots(ddf_sym, interval)
    if any(a < r < b for r in raices_ddf):
        raise ConvergenceError("f''(x) tiene raíces en el intervalo (punto de inflexión).")

    f0 = f.subs({x_var: x0}).n()
    ddf0 = ddf_sym.subs({x_var: x0}).n()
    if f0 * ddf0 <= 0:
        raise ConvergenceError("Condición de Fourier fallida: f(x0)*f''(x0) debe ser > 0 para garantizar que Newton no salte fuera.")

    return newton(f, x0, tol, itMax, minIt)

def newton(f: callable, x0: float, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[float, int]:
    nonNull(f, x0, tol, itMax, minIt)
    t = var('t')
    try:
        fSym = __fEval(f, t)
        dfSym = diff(fSym, t)
        

        g = lambda val: val - __fEval(f, val) / dfSym.subs(t=val)
    except Exception as e:
        raise ValueError(f"Error al derivar la función: {e}")
    root, it = fixedPointIteration(g, x0, tol, itMax, minIt)
    if isinstance(root, list):
        raise ValueError("Error raro. Contactar con soporte.")
    return root, it

def newtonSystem(F: callable, x0: vector, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[vector, int]:
    nonNull(F, x0, tol, itMax, minIt)
    n = len(x0)
    vars = [var(f'v{i}') for i in range(n)]

    fSym = F(vector(vars))
    JSym = jacobian(fSym, vars)
 

    def Nf(xVal):
        subs_dict = {vars[i]: RR(xVal[i]) for i in range(n)}
        try:
            FCurr = vector(RR, [expr.subs(subs_dict).n() for expr in fSym])

            JCurr = matrix(RR, n, n, [expr.subs(subs_dict).n() for expr in JSym.list()])
            delta = solvePLU(JCurr, -FCurr)
        except TypeError as e:
            raise ValueError(f"No se pudo evaluar la expresión a valor numérico: {e}")
        except ZeroDivisionError:
            raise ConvergenceError("Jacobiana singular detectada. Newton no convergerá.")
        except Exception as e:        
            print(f"DEBUG - Fallo en xVal: {xVal}")
            raise ValueError(f"Error en la evaluación numérica: {e}")
        return xVal + delta

    return fixedPointIteration(Nf, x0, tol, itMax, minIt)