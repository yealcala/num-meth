# num_meth/methods/nolinsys/newton.py

"""Provides some of the methods for newton / tangent method

- `newtonFourierSemilocal(f, interval, x0, tol, itMax, minIt)` - Applies newton iteration checking for Semilocal convergence (Fourier).
- `newton(f, x0, tol, itMax, minIt)` - Applies newton iteration.
- `newtonSystem(F, x0, tol, itMax, minIt)` - Applies newton iteration to a system.
"""
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
    """ Iterates over N(x_n) = x_{n+1} = x_n - f(x_n)/f'(x_n) checkinf first for Semilocal convergence (Fourier).

    Examples:
        >>> ec = x**6 - x - 1
        >>> show(newtonFourierSemilocal(ec, (1, 1.5), 1.3, tol=5e-5))
        

    Args:
        f (callable): any callable function.
        interval (tuple[float, float]): an interval.
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
        ConvergenceError: If convergence was not found for given interval or x0 point.
    """  
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
        raise ConvergenceError("Bolzano rule must work: f(a) and f(b) must have opposite signs.")


    raices_df = findGeneralRoots(dfSym, interval)
    if any(a < r < b for r in raices_df):
        raise ConvergenceError("f'(x) cannot have roots in the interval.")



    raices_ddf = findGeneralRoots(ddf_sym, interval)
    if any(a < r < b for r in raices_ddf):
        raise ConvergenceError("f''(x) cannot have roots in the interval.")

    f0 = f.subs({x_var: x0}).n()
    ddf0 = ddf_sym.subs({x_var: x0}).n()
    if f0 * ddf0 <= 0:
        raise ConvergenceError("f(x0)*f''(x0) must be positive.")

    return newton(f, x0, tol, itMax, minIt)

def newton(f: callable, x0: float, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[float, int]:
    """ Iterates over N(x_n) = x_{n+1} = x_n - f(x_n)/f'(x_n) 

    Examples:
        >>> f = e**x*log(x) + x**3 - 2
        >>> show(newton(f, 2, tol=1e-10))
        

    Args:
        f (callable): any callable function.
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
    nonNull(f, x0, tol, itMax, minIt)
    t = var('t')
    try:
        fSym = __fEval(f, t)
        dfSym = diff(fSym, t)
        

        g = lambda val: val - __fEval(f, val) / dfSym.subs(t=val)
    except Exception as e:
        raise ValueError(f"Couldn't differentiate funcion: {e}")
    root, it = fixedPointIteration(g, x0, tol, itMax, minIt)
    if isinstance(root, list):
        raise ValueError("Weird error. Please, consider contacting support.")
    return root, it

def newtonSystem(F: callable, x0: vector, tol: float = 1e-10, itMax: int = 100, minIt: int = 0) -> tuple[vector, int]:
    """ Iterates over N(x_n) = x_{n+1} = x_n - f(x_n)/f'(x_n) 

    Examples:
        >>> def F(v):
        >>>     return vector([v[0]**2 + v[1]**2 - 2*v[0] - 2*v[1] + 1, v[0] + v[1] - 2*v[0]*v[1] ])

        >>> show(newtonSystem(F, vector(RR, [2, 0])))

    Args:
        F (callable): a function as seen in the example.
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
        ConvergenceError: If divergence was found.
    """  
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
            raise ValueError(f"Cannot evaluate expression: {e}")
        except ZeroDivisionError:
            raise ConvergenceError("Jacobian is not inversible. It won't converge.")
        except Exception as e:        
            raise ValueError(f"Cannot evaluate expression: {e}")
        return xVal + delta

    return fixedPointIteration(Nf, x0, tol, itMax, minIt)